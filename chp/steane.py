"""
steane.py

"""

import sys
import os
import time
import json
from circuit import *
#from visualizer import Printer
import correction

class CircuitError(Exception): pass
class InvalidGateError(CircuitError): pass

class Code:

    """Defines constants and methods useful for dealing with the Steane Code"""

    code_name='Steane'

    """Gates that the code currently supports."""

    """Single qubit transversal gates."""
    single_qubit_unitaries = ["I","X","Y","Z","H","S"]
    """Single qubit non-clifford gates."""
    single_qubit_non_clifford_unitaries = ["T"]
    """Two qubit transversal gates."""
    two_qubit_unitaries = ["CX","XC","CZ","ZC","CS","SC"]
    """All transversal unitaries."""
    unitaries = single_qubit_unitaries + two_qubit_unitaries + single_qubit_non_clifford_unitaries
    """All preparations."""
    preparations =["PrepareZ", "PrepareX", "PrepareXPlus", "PrepareXMinus", 
                   "PrepareYPlus", "PrepareYMinus", "PrepareZPlus", "PrepareZMinus"]
    """All measurements."""
    measurements = ["MeasureX", "MeasureY", "MeasureZ",
                    "MeasureXDestroy", "MeasureYDestroy", "MeasureZDestroy"]

    """All error correcting routines."""
    error_corrections = ["EC_ShorCorrect","EC_SteaneCorrect","EC_KnillCorrect"]

    ecc_redundancy=3
    verify=False
    ancilla_parallel=True
    diVincenzo=True
    initial_I=False

    """Preparations and measurements have error correction built into them, so
    only gates to add error correction after are the unitaries."""
    error_corrected_gates = unitaries

    all_gates = unitaries + preparations + measurements + error_corrections

    stabilizer = [
            ["I","I","I","X","X","X","X"],
            ["I","X","X","I","I","X","X"],
            ["X","I","X","I","X","I","X"],
            ["I","I","I","Z","Z","Z","Z"],
            ["I","Z","Z","I","I","Z","Z"],
            ["Z","I","Z","I","Z","I","Z"]]

    stabilizer_alt = [
                        [('X',3), ('X',4), ('X',5), ('X',6)],
                        [('X',1), ('X',2), ('X',5), ('X',6)],
                        [('X',0), ('X',2), ('X',4), ('X',6)],
                        [('Z',3), ('Z',4), ('Z',5), ('Z',6)],
                        [('Z',1), ('Z',2), ('Z',5), ('Z',6)],
                        [('Z',0), ('Z',2), ('Z',4), ('Z',6)]]


    logical = {"X": ["X","X","X","I","I","I","I"],
           "Y": ["Y","Y","Y","I","I","I","I"],
               "Z": ["Z","Z","Z","I","I","I","I"]}

    Z_stabs = ['+XIXIXIX', '+IXXIIXX', '+IIIXXXX',
           '+ZZZIIII', '+IZZZZII', '+IIZIZZI', '+IIIZZZZ']
    Z_destabs = ['+ZIIIIII', '+IZIIIII', '+IIIZIII',
             '+IIXIXIX', '+IIIIXXI', '+IIIIIXX', '+IIIIIIX']
    X_stabs = ['+ZIZIZIZ', '+IZZIIZZ', '+IIIZZZZ',
           '+XXXIIII', '+IXXXXII', '+IIXIXXI', '+IIIXXXX']
    X_destabs = ['+XIIIIII', '+IXIIIII', '+IIIXIII',
             '+IIZIZIZ', '+IIIIZZI', '+IIIIIZZ', '+IIIIIIZ']

    stabs_dict = {'Z': (Z_stabs, Z_destabs),
                  'X': (X_stabs, X_destabs)}


    block_size = 7

    stabilizer_syndrome_dict = {
                    (0,0,0):["I","I","I","I","I","I","I"],
                    (0,0,1):["E","I","I","I","I","I","I"],
                    (0,1,0):["I","E","I","I","I","I","I"],
                    (0,1,1):["I","I","E","I","I","I","I"],
                    (1,0,0):["I","I","I","E","I","I","I"],
                    (1,0,1):["I","I","I","I","E","I","I"],
                    (1,1,0):["I","I","I","I","I","E","I"],
                    (1,1,1):["I","I","I","I","I","I","E"]}

    # look-up table in case the flag was triggered.
    # this applies when there is a single flag for the
    # 3 stabilizers.
    # syndromes (1,0,1) and (1,1,0) are only caused by
    # w-2 events.
    stabilizer_syndrome_dict_flag = {
                    (0,0,0):["I","I","I","I","I","I","I"],
                    (0,0,1):["I","E","E","I","I","I","I"],
                    (0,1,0):["E","I","E","I","I","I","I"],
                    (0,1,1):["I","I","E","I","I","I","I"],
                    (1,0,0):["I","I","I","E","I","I","I"],
                    (1,0,1):["I","I","I","I","E","I","I"],
                    (1,1,0):["I","I","I","I","I","E","I"],
                    (1,1,1):["I","I","I","I","I","I","E"]
                    }
                    
    total_lookup_table = {
                    0: stabilizer_syndrome_dict,
                    1: stabilizer_syndrome_dict_flag
                    }

    # Lookup table when we use 1 flag for each stabilizer
    # There is one dictionary for each one of the 4 possible
    # flag outcomes (no flags triggered, first, second, third)
    # The lookup tables are the same as the original with the
    # exception of the syndrome corresponding to a w-2 hook
    # error.  We assume that all w-2 error events will be
    # caused by a measurement error + data error.  Therefore,
    # for all other syndromes the correction is exactly the 
    # same as the the original one.
    syn_dict_flag0 = dict(stabilizer_syndrome_dict)
    syn_dict_flag0[(0,0,1)] = ["I","I","I","I","I","E","E"]
    syn_dict_flag1 = dict(stabilizer_syndrome_dict)
    syn_dict_flag1[(0,0,1)] = ["I","I","I","I","I","E","E"]
    syn_dict_flag2 = dict(stabilizer_syndrome_dict)
    syn_dict_flag2[(0,1,0)] = ["I","I","I","I","E","I","E"]

    total_lookup_table_one_flag = {
                                (0,0,0): stabilizer_syndrome_dict,
                                (1,0,0): syn_dict_flag0,
                                (0,1,0): syn_dict_flag1,
                                (0,0,1): syn_dict_flag2
                                }


    # stabilizers in CHP-friendly format
    
    stabilizer_CHP = [
                        '+XIXIXIX',
                        '+IXXIIXX',
                        '+IIIXXXX',
                        '+ZIZIZIZ',
                        '+IZZIIZZ',
                        '+IIIZZZZ'
                     ]

    destabilizer_CHP = [
                        '+ZIIIIII',
                        '+IZIIIII',
                        '+IIIZIII',
                        '+IIIIIXX',
                        '+IIIIXIX',
                        '+IIIIXXX'
                       ]

    stabilizer_CHP_Z = [
                        '+XIXIXIX',
                        '+IXXIIXX',
                        '+IIIXXXX',
                        '+ZZIZIIZ',
                        '+IZIIZIZ',
                        '+IIZZIIZ',
                        '+IIIZZZZ'
                       ]

    destabilizer_CHP_Z = [
                           '+ZIIIIII',
                           '+IZIIIII',
                           '+IIIZIII',
                           '+IIXIXIX',
                           '+IIIIXXI',
                           '+IIXIIII',
                           '+IIIIIXI'
                          ]

    stabilizer_CHP_X = [
                         '+XXXIIII',
                         '+IXXXXII',
                         '+IIXIXXI',
                         '+IIIXXXX',
                         '+ZIZIZIZ',
                         '+IZZIIZZ',
                         '+IIIZZZZ'
                       ]

    destabilizer_CHP_X = [
                           '+IIZIZIZ',
                           '+IIIIZZI',
                           '+IIIIIZZ',
                           '+IIIIIIZ',
                           '+XIIIIII',
                           '+IXIIIII',
                           '+IIIXIII'
                         ]

    
    stabilizer_logical_CHP = {'+Z': stabilizer_CHP_Z,
                              '+X': stabilizer_CHP_X,
                             }

    destabilizer_logical_CHP = {'+Z': destabilizer_CHP_Z,
                                '+X': destabilizer_CHP_X,
                               }




    @classmethod
    def decode_meas_Steane(cls, list_meas):
        '''
        returns the eigenvalues of the 3 stabilizers based on the measurement
        outcomes.
        list_meas is a list with 7 entries.
        s0 = c3 + c4 + c5 + c6
        s1 = c1 + c2 + c5 + c6
        s2 = c0 + c2 + c4 + c6
        '''
        s0 = (list_meas[3] + list_meas[4] + list_meas[5] + list_meas[6])%2
        s1 = (list_meas[1] + list_meas[2] + list_meas[5] + list_meas[6])%2
        s2 = (list_meas[0] + list_meas[2] + list_meas[4] + list_meas[6])%2

        return (s0, s1, s2)


    @classmethod
    def decode_syndrome_Steane_EC(cls, list_meas):
        '''
        '''
        s = cls.decode_meas_Steane(list_meas)
        return cls.stabilizer_syndrome_dict[s]



    @classmethod
    def parity_meas_Steane_EC(cls, list_meas):
        '''
        '''
        #print 'list meas =', list_meas
        b = sum(list_meas)%2
        s0, s1, s2 = cls.decode_meas_Steane(list_meas)
        if (s0 == 1) or (s1 == 1) or (s2 == 1):   
            return 1 - b
        else:                                     
            return b


    
    @classmethod
    def parity_check(cls, list_meas):
        '''
        '''
        s0, s1, s2 = cls.decode_meas_Steane(list_meas)
        if (s0 == 1) or (s1 == 1) or (s2 == 1):
            return 1
        else:
            return 0



    @classmethod
    def stabilizer_syndrome_to_error_circuits(cls, stab_type='X'):
        #given a syndrome of the form (+/-1,+/-1,+/-1), make a circuit with the correct error.

        corr_circ_list=[]
        for syndrome in stabilizer_syndrome_dict:
            corr_circ = Circuit()
            for index, err in enumerate(stabilizer_syndrome_dict[syndrome]):
                if err =='E':
                    corr_circ.add_gate_at([index],stab_type)
            corr_circ_list += [(syndrome,corr_circ)]
        return corr_circ_list


class Generator:
    """Steane code circuits
    """

    T_type='Nielsen and Chuang'

    @classmethod
    def create_encoded_circuit(cls, gate, parallel=True):
        """Creates an encoded circuit.  Valid circuits come from Code.
        """

        if type(gate) == str:
            gate_name = gate
        else:
            gate_name = gate.gate_name

        encoded_circ = None

        if (gate_name in Code.all_gates):
            if (gate_name in Code.single_qubit_unitaries):
                """Single qubit gates are transverse."""
                encoded_circ = Generator.single_qubit_gate(gate_name)
            elif (gate_name in Code.single_qubit_non_clifford_unitaries):
                """Single non-Clifford gate unitaries"""
                encoded_circ = Generator.non_clifford_single_gate(gate_name)
            elif (gate_name in Code.two_qubit_unitaries):
                """Two qubit gates are transverse."""
                encoded_circ = Generator.two_qubit_gate(gate_name)
            elif (gate_name=="EC_ShorCorrect"):
                """Cat based quantum error correction"""
                encoded_circ = Generator.shor_correct(parallel)
            elif (gate_name=="EC_SteaneCorrect"):
                """Steane based quantum error correction"""
                encoded_circ = Generator.steane_correct(parallel)
            elif (gate_name=="EC_KnillCorrect"):
                """Knill based quantum error correction"""
                encoded_circ = Generator.knill_correct()
            elif (gate_name[:7]=="Prepare"):
                """State preparations."""
                encoded_circ = Generator.pauli_prepare(gate_name[7])
            elif (gate_name[:7]=="Measure"):
                """Two types of measurements, destructive and non-desctructive."""
                if (gate_name[-7:]=="Destroy"):
                    encoded_circ = Generator.encoded_destructive_measurement(gate_name[7])
                else:
                    encoded_circ = Generator.encoded_pauli_measurement(gate_name[7])
            encoded_circ.parent_gate = gate
            return encoded_circ
        else:
            print(gate.gate_name)
            raise InvalidGateError

    @classmethod
    def single_qubit_gate(cls, gate_name=None):
        """Single qubit Steane code gates that are transversal.
        """
        circ = Circuit()
        for index in range(7):
            circ.add_gate_at([index],gate_name)
        return circ

    @classmethod
    def non_clifford_single_gate(cls, gate_name=None):

        if(gate_name=="T"):
            if cls.T_type == 'Nielsen and Chuang' :
                """
                From Nielsen and Chuang p485.
                The following is equivalent to a fault tolerant T gate:

                |phi>------------------------------------------O--_ _-------|SX|-----T|phi>
                                                               |  _X_        ||
                       |0>----|H|----|T|---|exp(-i pi/4)SX|----.--  -----|Measure|
                                                 |
                                         |+>-----.----|MeasureX|
                """
                redundancy = 1

                magic_state_circ =  cls._magic_state_prep_and_verify_Nielsen_and_Chuang(redundancy)

                measure_circ = Circuit()

                measure_circ.join_circuit_at(list(range(7,14)),magic_state_circ)

                for i in range(7):
                    measure_circ.add_gate_at([i+7,i],'CX')
                    measure_circ.add_gate_at([i],"MeasureZDestroy")

                switching_qubits = measure_circ._get_qubits(list(range(14)))
                permute_gate = Permute_Gate('permute',switching_qubits,switching_qubits[7:14]+switching_qubits[:7])

                correct_T_syndrome = Circuit()
                for i in range(7):
                    correct_T_syndrome.add_gate_at([7+i],"X")
                    correct_T_syndrome.add_gate_at([7+i],"S")
                correct_T_syndrome_gate = Correction_Gate(gate_name= "Tcorrect", circuit_list=[correct_T_syndrome])

                #Put it all together:

                T_circuit = Circuit()
                T_circuit.join_circuit(measure_circ)

                T_circuit.add_gate_at(list(range(14)),permute_gate)

                T_circuit.to_ancilla(list(range(7,14)))

                #Printer.circuit_console_drawing(T_circuit)

                return T_circuit

            elif T_Type== 'Aliferis Magic State':
                pass

            elif T_Type== 'ORAQL GFI':
                pass


        else:
            raise InvalidGateError

    @classmethod
    def _magic_state_prep_and_verify_Nielsen_and_Chuang(cls, redundancy=1):

        #create magic state:
        magic_prep_circ = Circuit()

        magic_prep_circ.join_circuit(Generator.FT_encoded_zero())

        for i in range(7):
            magic_prep_circ.add_gate_at([i],"H")

        for i in range(7):
            magic_prep_circ.add_gate_at([i],"T")

        #measure magic state
        controls = []
        for i_redund in range(redundancy):

            synd_circ= Circuit()

            enc_plus = Generator.create_encoded_plus()

            synd_circ.join_circuit_at(list(range(7,14)),enc_plus)

            for i in range(7):
                synd_circ.add_gate_at([i,7+i],"CX")
                synd_circ.add_gate_at([i,7+i],"CS")
                synd_circ.add_gate_at([7+i],"T")
                synd_circ.add_gate_at([7+i],"S")
                synd_circ.add_gate_at([7+i],"Z")

            #(The operator ZST is T^dag and is meant to do the controlled-e^(-i pi/4) gate)

            for i in range(7):
                synd_circ.add_gate_at([7+i],'X')

            controls += [Classical_Control(function= 'parity', answer = 1,
                        syndrome_qubits= synd_circ._get_qubits(list(range(7,14))))]

            synd_circ.to_ancilla(list(range(7,14)))

            magic_prep_circ.join_circuit(synd_circ)

        control = Classical_Control(function = 'redundancy', answer = 1, controls= controls)

        magic_gate = Verify_Gate(gate_name='MagicStateTheta', circuit_list=[magic_prep_circ])

        #return magic_gate.circuit_wrap()
        return magic_prep_circ

    @classmethod
    def two_qubit_gate(cls, gate=None):
        """Two qubit Steane code gates that are transversal.
        """
        circ = Circuit()
        for index in range(7):
            circ.add_gate_at([index,index+7],gate)
        return circ

    @classmethod
    def encoded_destructive_measurement(cls, gate_name=None, meas_errors=True):
        """This defines a destructive measurement circuit of an encoded pauli gate
        of a CSS code.
        Note that in order to make this fault-tolerant one must use this carefully.
        In particulary one should not rely on the product of the measurement gates.
        Instead one should run a simple redundancy code check on these outputs and
        then interpret the measurement as the product of the measurment gates. """

        if not (gate_name=='X' or gate_name == 'Z'):
            print('Can only encode X or Z destructive measurements this way.')
            return False

        n_data = Code.block_size

        I_error = 'Im' + gate_name
        m_circ=Circuit()
        for index in range(n_data):
            if meas_errors:
                m_circ.add_gate_at([index], I_error)
            m_circ.add_gate_at([index],'Measure'+gate_name+'Destroy')

        return m_circ


    # At the moment this is not implemented. When it is, the decoder will have to be
    # modified. -Lukas (12-4-2011)
    @classmethod
    def encoded_pauli_measurement(cls, gate_name=None, redundancy=1, return_gate = False):
        '''
        non-destructively measurement of the logical qubit using a 3-cat state.
        '''
        correct = correction.Cat_Correct.create_cat_3

        n_data = Code.block_size
        n_ancilla= len(correct().qubits())-n_data

        circ = Circuit()

        logical_gate=Code().logical[gate_name]

        for i_redund in range(redundancy-1):
            for index in range(n_data):

                anc_start=n_data+i_redund*n_ancilla
                anc_end=n_data+(i_redund+1)*n_ancilla

                corr_circ = correct()
                circ.join_circuit_at(list(range(anc_start,anc_end)),corr)

                if(logical_gate==gate_name):
                    circ.add_gate_at([anc_start,index],'C'+gate_name)
                    circ.add_gate_at([index],'Measure'+gate_name+'Destroy')

        gate = Encoded_Gate(gate_name='Encoded Measure'+gate_name,circuit_list=[circ])

        if return_gate==False:  return gate.circuit_wrap()
        if return_gate==True:   return gate


    @classmethod
    def pauli_prepare(cls, basis):
        if basis == 'Z':
            return correction.Steane_Correct.FT_encoded_zero_Steane()[0]
        elif basis == 'X':
            return correction.Steane_Correct.FT_encoded_plus_Steane()[0]





    # Error correction methods

    @classmethod
    def shor_correct(cls, parallel):
        redundancy = Code.ecc_redundancy
        verify = Code.verify
        ancilla_parallel = parallel
        diVincenzo = Code.diVincenzo
        initial_I = Code.initial_I
        return correction.Cat_Correct.cat_syndrome_4(Code.stabilizer,redundancy,
                                                     verify, ancilla_parallel,
                                                     diVincenzo, initial_I)

    @classmethod
    def steane_correct(cls, parallel):
        #redundancy=Code.ecc_redundancy
        # For Steane EC, we don't need to repeat 3 times to make it FT.
        #redundancy = 1
        #steane_module=sys.modules[__name__]
        #return correction.Steane_Correct.steane_syndrome(ecc=steane_module,
        #                                                 redundancy=redundancy,
        ancilla_parallel = parallel
        return correction.Steane_Correct.steane_syndrome(ancilla_parallel)

    @classmethod
    def knill_correct(cls):
        steane_module=sys.modules[__name__]
        return correction.Knill_Correct.knill_syndrome(ecc=steane_module)



    @classmethod
    def physical_CNOT_ion_trap(cls, CNOT_index=0):
        '''
        Circuit for a fully scheduled transversal CNOT on an ion trap based
        on Alejandro Bermudez's work and eQual's experimental parameters
        
        The basic idea is that 2 logical qubits involved live on different
        arms of a Y-junction trap.  We cannot perform entangling MS gates
        between non-nearest neighbors, so we have to shuttle ions back and
        forth.

        Duration assumptions (in arbitrary units; only ratios matter):
        Shuttle = 1; splitting/merging = 3; rotation = 2; 1q gate = 0;
        2q gate = 2; multi-qubit gate = 2; measurement = 3; prep = 1;
        cooling = 10.
        The duration of the transport through the junction is an 
        independent parameter.
        '''

        if CNOT_index == 0:
            waiting_times = [6,5,1]
            #n_cross = 3
            n_cross = 9
        elif CNOT_index == 4:
            waiting_times = [8,5,1]
            #n_cross = 4
            n_cross = 14
        else:
            waiting_times = [5,5,1]
            n_cross = 0

        CNOT_circ = Circuit()
        
        # Initial shuttling of group of 4 ions
        for q in range(7):
            for t in range(waiting_times[0]):
                CNOT_circ.add_gate_at([q], 'I_idle')
                CNOT_circ.add_gate_at([q+7], 'I_idle')
            for t in range(n_cross):
                CNOT_circ.add_gate_at([q], 'I_cross')
                CNOT_circ.add_gate_at([q+7], 'I_cross')

        # Y rotation on ctrl qubit
        CNOT_circ.add_gate_at([CNOT_index], 'RY +')

        # Shuttling and cooling before MS gate
        for q in range(7):
            for t in range(waiting_times[1]):
                CNOT_circ.add_gate_at([q], 'I_idle')
                CNOT_circ.add_gate_at([q+7], 'I_idle')

        # MS gate and X rotations
        # Duration of MS gate = 1 time unit
        active_qubits = [CNOT_index, CNOT_index+7]
        CNOT_circ.add_gate_at(active_qubits, 'MS')
        idle_qubits = [i for i in range(14) if i not in active_qubits]
        for i in idle_qubits:
            CNOT_circ.add_gate_at([i], 'I_idle')
        
        CNOT_circ.add_gate_at([CNOT_index], 'RX -')
        CNOT_circ.add_gate_at([CNOT_index+7], 'RX -')

        # Shuttling before second Y rotation
        for i in range(7):
            for t in range(waiting_times[2]):
                CNOT_circ.add_gate_at([i], 'I_idle')
                CNOT_circ.add_gate_at([i+7], 'I_idle')

        # Second Y rotation on ctrl qubit
        CNOT_circ.add_gate_at([CNOT_index], 'RY -')
        
        return CNOT_circ


    
    @classmethod
    def transversal_CNOT_ion_trap(cls, encoded_CNOTs=True, encoded_total_circ=False):
        '''
        All the 7 physical CNOTs
        '''

        total_circ = Circuit()
        for CNOT_index in range(7):
            CNOT_circ = Generator.physical_CNOT_ion_trap(CNOT_index)
            if encoded_CNOTs:
                CNOT_name = 'CNOT_%i' %CNOT_index
                CNOT_circ = Encoded_Gate(CNOT_name, [CNOT_circ]).circuit_wrap()
            total_circ.join_circuit(CNOT_circ)

        # Add the final re-ordering
        final_circ = Circuit()
        for q in range(14):
            for t in range(6):
                final_circ.add_gate_at([q], 'I_idle')
            #n_cross_final = 3
            n_cross_final = 9
            for t in range(n_cross_final):
                final_circ.add_gate_at([q], 'I_cross')

        if encoded_CNOTs:
            final_name = 'Final_reordering'
            final_circ = Encoded_Gate(final_name, [final_circ]).circuit_wrap()
        total_circ.join_circuit(final_circ)

        if encoded_total_circ:
            total_circ = Encoded_Gate('Transversal_CNOT', [total_circ]).circuit_wrap()

        return total_circ


    # The following methods were moved to correction.Steane:

    #@classmethod
    #def create_encoded_zero(cls):
    #   """Creates an encoded zero for explicit use only in Steane."""
    #   enc_zero = Circuit()
    #   for index in [0,1,3]:
    #       enc_zero.add_gate_at([index],"PrepareXPlus")
    #   for index in [2,4,5,6]:
    #       enc_zero.add_gate_at([index],"PrepareZPlus")
    #   cnot_sequence=[[0,2],[3,5],[1,6],[0,4],[3,6],[1,5],[0,6],[1,2],[3,4]]
    #   for indexes in cnot_sequence:
    #       enc_zero.add_gate_at(indexes,"CX")

    #   return enc_zero


    #@classmethod
    #def create_encoded_plus(cls):
    #   """Creates an encoded plus for explicit use only in Steane code.
    #   Almost identical to create_econded_zero; the only difference being CNOTs
    #   are flipped, 0's are +'s, and +'s are 0's."""
    #   enc_plus = Circuit()
    #   for index in [0,1,3]:
    #       enc_plus.add_gate_at([index],"PrepareZPlus")
    #   for index in [2,4,5,6]:
    #       enc_plus.add_gate_at([index],"PrepareXPlus")
    #   cnot_sequence=[[2,0],[5,3],[6,1],[4,0],[6,3],[5,1],[6,0],[2,1],[4,3]]
    #   for indexes in cnot_sequence:
    #       enc_plus.add_gate_at(indexes,"CX")

    #   enc_plus2=Circuit()
    #   for i in range(7):
    #       enc_plus2.add_gate_at([i],'Z')

    #   return enc_plus


    #@classmethod
    #def FT_encoded_zero(cls):

    #   purify_circ = cls.create_encoded_zero()
    #   enc_zero_anc= cls.create_encoded_zero()
    #   purify_circ.join_circuit_at(range(7,14), enc_zero_anc)

    #   for i in range(7):
    #       purify_circ.add_gate_at([i,7+i],'CX')

        #purify_circ._get_qubits(range(7))
        #measure and correct:

    #   for i in range(7):
    #       purify_circ.add_gate_at([7+i],'MeasureZDestroy')

        #control = Classical_Control(function='parity',
        #       syndrome_qubits = purify_circ._get_qubits(range(7,14)),
        #       target_qubits   = purify_circ._get_qubits(range(7)))

        #corr_circ = Circuit()
        #for i in range(7):
        #   corr_circ.add_gate_at([i],'X')

        #gate = Correction_Gate('Correct_To_Logical_Zero', control, [corr_circ])

        #purify_circ.add_gate(gate)

    #   purify_circ.to_ancilla(range(7,14))

    #   return purify_circ


    #@classmethod
    #def FT_encoded_plus(cls):
    #   purify_circ = cls.create_encoded_plus()
    #   enc_plus_anc= cls.create_encoded_plus()
    #   purify_circ.join_circuit_at(range(7,14), enc_plus_anc)

    #   for i in range(7):
    #       purify_circ.add_gate_at([7+i,i],'CX')

        #measure and correct:

    #   for i in range(7):
    #       purify_circ.add_gate_at([7+i],'MeasureXDestroy')

        #control = Classical_Control(function='parity',
        #       syndrome_qubits = purify_circ._get_qubits(range(7,14)),
        #       target_qubits   = purify_circ._get_qubits(range(7)))

        #corr_circ = Circuit()
        #for i in range(7):
        #   corr_circ.add_gate_at([i],'Z')

        #gate = Correction_Gate('Correct_To_Logical_Zero', control, [corr_circ])

        #purify_circ.add_gate(gate)

        #purify_circ.to_ancilla(range(7,14))

    #   return purify_circ
