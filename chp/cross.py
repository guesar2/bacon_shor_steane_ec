'''
cross.py
'''

import sys
import os
import time
import json
from circuit import *
import correction


class Code:
    '''
    Defines constants and methods useful for dealing with the Cross 7-qubit code
    '''

    code_name = 'Cross'

    stabilizer = [
                   ['X','I','I','I','X','I','I'],
                   ['I','X','I','I','X','I','I'],
                   ['I','I','X','I','I','X','I'],
                   ['I','I','I','X','I','I','X'],
                   ['I','I','Z','Z','I','Y','Y'],
                   ['Z','Z','Z','X','Z','Z','I']]
    
    stabilizer_alt = [
                        [('X',0), ('X',4)],
                        [('X',1), ('X',4)],
                        [('X',2), ('X',5)],
                        [('X',3), ('X',6)],
                        [('Z',2), ('Z',3), ('Y',5), ('Y',6)],
                        [('Z',0), ('Z',2), ('X',3), ('Z',1), ('Z',4), ('Z',5)]]

    # stabilizers in CHP-friendly format
    stabilizer_CHP = [
                      '+XIIIXII',
                      '+IXIIXII',
                      '+IIXIIXI',
                      '+IIIXIIX',
                      '+IIZZIYY',
                      '+ZZZIZZX']
    
    destabilizer_CHP = [
                        '+ZIIIIII',
                        '+IZIIIII',
                        '+IIZIIII',
                        '+IIIZIII',
                        '+IIZIIZI',
                        '+IIZZIZZ']
    
    stabilizer_logical_CHP = {'+Z': '+ZZIIZII',
                              '-Z': '-ZZIIZII'}  

    destabilizer_logical_CHP = {'+Z': '+IIZZXZZ',
                                '-Z': '+IIZZXZZ'}

    # The basic lookup table is the one presented in Andrew
    # Cross's notes.  It includes the 21 single-qubit Pauli
    # errors plus 4 higher-weight errors that correspond
    # to single-qubit errors on the ancilla propagating
    # to the data.
    basic_lookup_table = {
                            0:  ['I','I','I','I','I','I','I'],
                            32: ['Z','I','I','I','I','I','I'],
                            16: ['I','Z','I','I','I','I','I'],
                            8:  ['I','I','Z','I','I','I','I'],
                            5:  ['I','I','I','Z','I','I','I'],
                            48: ['I','I','I','I','Z','I','I'],
                            10: ['I','I','I','I','I','Z','I'],
                            6:  ['I','I','I','I','I','I','Z'],
                            1:  ['X','I','I','I','I','I','I'],
                            3:  ['I','I','X','I','I','I','I'],
                            2:  ['I','I','I','X','I','I','I'],
                            33: ['Y','I','I','I','I','I','I'],
                            17: ['I','Y','I','I','I','I','I'],
                            11: ['I','I','Y','I','I','I','I'],
                            7:  ['I','I','I','Y','I','I','I'],
                            49: ['I','I','I','I','Y','I','I'],
                            9:  ['I','I','I','I','I','Y','I'],
                            4:  ['I','I','I','I','I','I','Y'],
                            13: ['I','I','Z','Z','I','I','I'],
                            40: ['Z','I','Z','I','I','I','I'],
                            42: ['Z','I','Z','X','I','I','I'],
                            58: ['I','I','I','I','Z','Z','I']
                           }


    # The complete lookup table is imported from the json file.
    # After the import, we convert each key from a string to
    # an integer.
    # The json file was generated with the script 'Cross_decoding.py'.
    # The idea is to add all the errors in the basic lookup table
    # until we have obtained the 64 possible syndromes. 
    table_json_filename = "complete_lookup_table.json"
    table_json_path = os.path.join(os.path.dirname(__file__), table_json_filename)
    table_json = open(table_json_path, 'r')
    complete_lookup_table_str = json.load(table_json)
    table_json.close()
    complete_lookup_table = {}
    for key in complete_lookup_table_str:
        complete_lookup_table[int(key)] = complete_lookup_table_str[key]


