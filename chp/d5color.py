import sys
import os
import json
from .circuit import *
from . import correction


class Code:

    code_name = 'd5color'

    # The weight-8 stabilizer always comes first by convention
    stabilizer = [
    ['I', 'I', 'X', 'X', 'I', 'X', 'X', 'I', 'I', 'X', 'X', 'I', 'I', 'X', 'X', 'I', 'I'],
    ['X', 'X', 'X', 'X', 'I', 'I', 'I', 'I', 'I', 'I', 'I', 'I', 'I', 'I', 'I', 'I', 'I'],
    ['X', 'I', 'X', 'I', 'X', 'X', 'I', 'I', 'I', 'I', 'I', 'I', 'I', 'I', 'I', 'I', 'I'],
    ['I', 'I', 'I', 'I', 'X', 'X', 'I', 'I', 'X', 'X', 'I', 'I', 'I', 'I', 'I', 'I', 'I'],
    ['I', 'I', 'I', 'I', 'I', 'I', 'I', 'I', 'X', 'X', 'I', 'I', 'X', 'X', 'I', 'I', 'I'],
    ['I', 'I', 'I', 'I', 'I', 'I', 'X', 'X', 'I', 'I', 'X', 'X', 'I', 'I', 'I', 'I', 'I'],
    ['I', 'I', 'I', 'I', 'I', 'I', 'I', 'I', 'I', 'I', 'X', 'X', 'I', 'I', 'X', 'X', 'I'],
    ['I', 'I', 'I', 'I', 'I', 'I', 'I', 'X', 'I', 'I', 'I', 'X', 'I', 'I', 'I', 'X', 'X'],
    ['I', 'I', 'Z', 'Z', 'I', 'Z', 'Z', 'I', 'I', 'Z', 'Z', 'I', 'I', 'Z', 'Z', 'I', 'I'],
    ['Z', 'Z', 'Z', 'Z', 'I', 'I', 'I', 'I', 'I', 'I', 'I', 'I', 'I', 'I', 'I', 'I', 'I'],
    ['Z', 'I', 'Z', 'I', 'Z', 'Z', 'I', 'I', 'I', 'I', 'I', 'I', 'I', 'I', 'I', 'I', 'I'],
    ['I', 'I', 'I', 'I', 'Z', 'Z', 'I', 'I', 'Z', 'Z', 'I', 'I', 'I', 'I', 'I', 'I', 'I'],
    ['I', 'I', 'I', 'I', 'I', 'I', 'I', 'I', 'Z', 'Z', 'I', 'I', 'Z', 'Z', 'I', 'I', 'I'],
    ['I', 'I', 'I', 'I', 'I', 'I', 'Z', 'Z', 'I', 'I', 'Z', 'Z', 'I', 'I', 'I', 'I', 'I'],
    ['I', 'I', 'I', 'I', 'I', 'I', 'I', 'I', 'I', 'I', 'Z', 'Z', 'I', 'I', 'Z', 'Z', 'I'],
    ['I', 'I', 'I', 'I', 'I', 'I', 'I', 'Z', 'I', 'I', 'I', 'Z', 'I', 'I', 'I', 'Z', 'Z']
        ]


    # The lookup table is imported from the json file.
    table_json_filename = 'complete_lookup_table_d5.json'
    table_json = open(table_json_filename, 'r')
    lookuptable_str = json.load(table_json)
    table_json.close()


    stabilizer_CHP_X = [
                        '+XXXXIIIIIIIIIIIII',
                        '+IXIXXXIIIIIIIIIII',
                        '+IIXXXIIIIXIIIXIII',
                        '+IIIIXXIIXXIIIIIII',
                        '+IIIIIIXXIIXXIIIII',
                        '+IIIIIIIXXXIXIIXII',
                        '+IIIIIIIIXXIIXXIII',
                        '+IIIIIIIIIIXXIIXXI',
                        '+IIIIIIIIIIIIXXXXX',
                        '+ZZZZIIIIIIIIIIIII',
                        '+IZIZZZIIIIIIIIIII',
                        '+IIZZIZZIIZZIIZZII',
                        '+IIIIZZIIZZIIIIIII',
                        '+IIIIIIZZIIZZIIIII',
                        '+IIIIIIIZIIIZIIIZZ',
                        '+IIIIIIIIZZIIZZIII',
                        '+IIIIIIIIIIZZIIZZI'
                        ]

    destabilizer_CHP_X = [
                            '+IIZIIIIIIIIIZZIII',
                            '+IIZZIIZIZZZIZZIZZ',
                            '+IIIIIIIIIIIIZZIII',
                            '+IIIIIIIIZIIIZIZZZ',
                            '+IIIIIIIIIIIZIIZIZ',
                            '+IIIIIIIIIIIIIIZZI',
                            '+IIIIIIIIIIIIZIIIZ',
                            '+IIIIIIIIIIIIIIIZZ',
                            '+IIIIIIIIIIIIIIIIZ',
                            '+XIIIIIIIIIIIIIIII',
                            '+XXIIIIIIIIIIIIIII',
                            '+IIIIXXIIIIIIIIIII',
                            '+XXIIXIIIIIIIIIIII',
                            '+IXIXIIXIIIIIIIIII',
                            '+IXIXIIXXIIIIIIIII',
                            '+XIIXXIIIIXIIIIIII',
                            '+IIIIIIXIIIXIIIIII'
                          ]


    stabilizer_alt = [
                        [('X',2),('X',3),('X',6),('X',5),('X',9),('X',10),('X',14),('X',13)],
                        [('X',0),('X',1),('X',2),('X',3)],
                        [('X',0),('X',2),('X',4),('X',5)],
                        [('X',4),('X',5),('X',8),('X',9)],
                        [('X',8),('X',9),('X',12),('X',13)],
                        [('X',6),('X',7),('X',10),('X',11)],
                        [('X',10),('X',11),('X',14),('X',15)],
                        [('X',7),('X',11),('X',15),('X',16)],
                        [('Z',2),('Z',3),('Z',6),('Z',5),('Z',9),('Z',10),('Z',14),('Z',13)],
                        [('Z',0),('Z',1),('Z',2),('Z',3)],
                        [('Z',0),('Z',2),('Z',4),('Z',5)],
                        [('Z',4),('Z',5),('Z',8),('Z',9)],
                        [('Z',8),('Z',9),('Z',12),('Z',13)],
                        [('Z',6),('Z',7),('Z',10),('Z',11)],
                        [('Z',10),('Z',11),('Z',14),('Z',15)],
                        [('Z',7),('Z',11),('Z',15),('Z',16)]
                     ]


    lookup_folder = './color_code_with_flags/lookup_tables_color_d5/'
    lookup_filenames = os.listdir(lookup_folder)
    all_lookups = {}
    for lookup_filename in lookup_filenames:
        total_filename = lookup_folder + lookup_filename
        lookup_file = open(total_filename, 'r')
        lookup_dict = json.load(lookup_file)
        lookup_file.close()
        all_lookups[lookup_filename[:-5]] = lookup_dict
