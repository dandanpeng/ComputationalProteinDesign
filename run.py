#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 09:12:39 2019

@author: pengdandan
"""
from function import *
from seqpred import *

ids = ['1aba', '1bxv', '1by2',  '1hyp', '1opc', 
       '1tmy', '2acy', '2mcm', '3e0e', '3k63']

score = 0
for i in ids: 
    protein = Protein(i)
    match = Matches(protein.match)
    genetic_algorithm(protein, match, 100, 50, 2, 0.03, 2)
    

