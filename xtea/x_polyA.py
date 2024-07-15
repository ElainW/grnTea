##11/01/2018
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu

# YW notes: defined a python class: PolyA

####canonical signals: AATAAA / TTTATT, ATTAAA/TTTAAT
####variations: AGTAAA, TATAAA, CATAAA, GATAAA, AATATA, AATACA, AATAGA, AAAAAG, ACTAAA, AAGAAA, AATGAA, TTTAAA, AAAACA,
# GGGGCT
import re

class PolyA():


    #check whether contain enough A or T
    # YW 2020/08/09 added min_ratio
    def contain_enough_A_T(self, seq, seq_len, min_ratio):
        n_A=0
        n_T=0
        for ch in seq:
            if ch == 'A' or ch == 'a':
                n_A+=1
            if ch == 'T' or ch == 't':
                n_T+=1

        if n_A>=float(min_ratio*seq_len) or n_T>=float(min_ratio*seq_len):
            return True
        return False
    
