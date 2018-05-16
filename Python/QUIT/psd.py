#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Implements a QUIT sequences for easy prototyping of sequences in Python
"""

import json
import os

class SequenceGroup(object):
    """
    Implementation of a Sequence Group for composite sequences such as mcDESPOT
    """

    def __init__(self, sequences=None):
        self.seq_group = {
            'SequenceGroup':{
                'sequences':[]
                }
            }

        if len(sequences) > 0:
            for seq in sequences:
                self.add_sequence(seq)

    def add_sequence(self, seq):
        d = seq.seq
        
        # Safety check (Should never occur)
        if len(d.keys()) > 1:
            print('Sequence can only have one key')
        
        # if list(d.keys())[0] in self.seq_group.keys():
        #     print('Sequence is already in the group. Not added')
        # else:
        
        self.seq_group['SequenceGroup']['sequences'].append(d)

    def save_json(self, filename):
        if os.path.exists(filename):
            raise FileExistsError('The file already exist')
        else:
            with open(filename, 'w') as fp:
                json.dump(self.seq_group, fp)

class Sequence(object):

    def __init__(self, name):
        self.name = name
        self.params = {}
        self.seq = {name:self.params}
        
    def save_json(self, filename):
        if os.path.exists(filename):
            raise FileExistsError('The file already exist')
        else:
            with open(filename, 'w') as fp:
                json.dump(self.seq, fp)

    def __repr__(self):
        return "<{} Sequence Object: {}>".format(self.name, self.seq)

    def __str__(self):
        return "<{} Sequence Object: {}>".format(self.name, self.seq)

class AFI(Sequence):

    def __init__(self, TR1, TR2, FA):
        super().__init__('AFI')
        self.params['TR1'] = TR1
        self.params['TR2'] = TR2
        self.params['FA'] = FA

class CASL(Sequence):

    def __init__(self, TR, label_time, post_label_delay):
        super().__init__('CASL')
        self.params['TR'] = TR
        self.params['label_time'] = label_time
        self.params['post_label_delay'] = post_label_delay

class MPRAGE(Sequence):

    def __init__(self, TR, TI, TD, eta, FA, etl, k0):
        super().__init__('MPRAGE')
        self.params['TR'] = TR
        self.params['TI'] = TI
        self.params['TD'] = TD
        self.params['eta'] = eta
        self.params['FA'] = FA
        self.params['etl'] = etl
        self.params['k0'] = k0

class MTSat(Sequence):

    def __init__(self, FA, TR, sat_f0, sat_angle,):
        super().__init__('MTSat')
        self.params['FA'] = FA
        self.params['TR'] = TR
        self.params['sat_f0'] = sat_f0
        self.params['sat_angle'] = sat_angle

class MultiEcho(Sequence):
    
    def __init__(self, TR, TE1, ESP, ETL):
        super().__init__('MultiEcho')
        self.params['TR'] = TR
        self.params['TE1'] = TE1
        self.params['ESP'] = ESP
        self.params['ETL'] = ETL

class SPGR(Sequence):
    
    def __init__(self, FA, TR):
        super().__init__('SPGR')
        self.params['FA'] = FA
        self.params['TR'] = TR

class SSFP(Sequence):
    
    def __init__(self, TR, FA, phaseinc):
        super().__init__('SSFP')
        self.params['TR'] = TR
        self.params['FA'] = FA
        self.phaseinc = phaseinc