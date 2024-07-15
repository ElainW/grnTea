##11/22/2017
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu

####Given rmsk annotation file, load the file into memory
####Check whether a given position is within a repeat region
####Parse out the given seq regions in the reference

####Revision on 03/04/2019
####Use interval tree to replace binary search

# YW notes: define the following class: XAnnotation

import os
import sys
import pysam
import global_values
from cmd_runner import *
from intervaltree import IntervalTree

####
class XAnnotation():
    def __init__(self, sf_rmsk):
        self.sf_annotation = sf_rmsk
        ##in format: m_rmsk_annotation[chrm][start_pos].append((end_pos, b_rc, sub_type, csn_start, csn_end))
        self.m_rmsk_annotation = {}
        self.m_rmsk_chrm_regions = {}
        self.b_with_chr = True
        self.m_interval_tree = {}  # by chrm

    # This indicates whether the chromosomes is in format "chr1" or "1"
    def set_with_chr(self, b_with_chr):
        self.b_with_chr = b_with_chr

    ## "self.b_with_chr" is the format gotten from the alignment file
    ## all other format should be changed to consistent with the "self.b_with_chr"
    def _process_chrm_name(self, chrm):
        b_chrm_with_chr = False
        if len(chrm) > 3 and chrm[:3] == "chr":  ##Here remove the "chr"
            b_chrm_with_chr = True

        # print chrm, self.b_with_chr, b_chrm_with_chr #################################################################

        if self.b_with_chr == True and b_chrm_with_chr == True:
            return chrm
        elif self.b_with_chr == True and b_chrm_with_chr == False:
            return "chr" + chrm
        elif self.b_with_chr == False and b_chrm_with_chr == True:
            return chrm[3:]
        else:
            return chrm


    '''
    Each record of the annotation file is in the format:
        1306    = Smith-Waterman score of the match, usually complexity adjusted
            The SW scores are not always directly comparable. Sometimes
            the complexity adjustment has been turned off, and a variety of
            scoring-matrices are used.
        15.6    = % substitutions in matching region compared to the consensus
        6.2     = % of bases opposite a gap in the query sequence (deleted bp)
        0.0     = % of bases opposite a gap in the repeat consensus (inserted bp)
        HSU08988 = name of query sequence
        6563    = starting position of match in query sequence
        7714    = ending position of match in query sequence
        (22462) = no. of bases in query sequence past the ending position of match
        C       = match is with the Complement of the consensus sequence in the database
        MER7A   = name of the matching interspersed repeat
        DNA/MER2_type = the class of the repeat, in this case a DNA transposon 
                fossil of the MER2 group (see below for list and references)
        (0)     = no. of bases in (complement of) the repeat consensus sequence 
                prior to beginning of the match (so 0 means that the match extended 
                all the way to the end of the repeat consensus sequence)
        2418    = starting position of match in database sequence (using top-strand numbering)
        1465    = ending position of match in database sequence
    '''
####
    # 7288   1.3  0.1  0.1  chr21     32213340 32214177 (15915718) +  L1HS LINE/L1 5312 6155    (0) 4482983
    # 4741   0.6  0.0  0.0  chr21     17916705 17917238 (30212657) C  L1HS LINE/L1 (0) 6155   5622 4462908
    def load_rmsk_annotation_no_extnd(self):
        with open(self.sf_annotation) as fin_rmsk:
            for line in fin_rmsk:
                fields = line.split()
                tmp_chrm = fields[4]
                chrm = self._process_chrm_name(tmp_chrm)
                start_pos = int(fields[5])
                end_pos = int(fields[6])
                b_rc = False
                if fields[8] == "C":
                    b_rc = True
                sub_type = fields[9]
                csn_start = -1
                csn_end = -1
                if b_rc == True:
                    csn_start = int(fields[13])
                    csn_end = int(fields[12])
                else:
                    csn_start = int(fields[11])
                    csn_end = int(fields[12])

                if chrm not in self.m_rmsk_annotation:
                    self.m_rmsk_annotation[chrm] = {}
                if start_pos in self.m_rmsk_annotation[chrm]:  ##this is not allowed!!!!
                    print(("Position {0}:{1} has more than  1 annotation!!!!".format(chrm, start_pos)))

                extd_start_pos = start_pos
                if extd_start_pos not in self.m_rmsk_annotation[chrm]:
                    self.m_rmsk_annotation[chrm][extd_start_pos] = []
                self.m_rmsk_annotation[chrm][extd_start_pos].append(
                    (end_pos, b_rc, sub_type, csn_start, csn_end))
                

    # 7288   1.3  0.1  0.1  chr21     32213340 32214177 (15915718) +  L1HS LINE/L1 5312 6155    (0) 4482983
    # 4741   0.6  0.0  0.0  chr21     17916705 17917238 (30212657) C  L1HS LINE/L1 (0) 6155   5622 4462908
    def load_rmsk_annotation_with_extnd_with_lenth_cutoff(self, i_extnd, i_min_len):
        with open(self.sf_annotation) as fin_rmsk:
            for line in fin_rmsk:
                fields = line.split()
                tmp_chrm = fields[4]
                chrm = self._process_chrm_name(tmp_chrm)
                start_pos = int(fields[5])
                end_pos = int(fields[6])
                b_rc = False
                if fields[8] == "C":
                    b_rc = True
                sub_type = fields[9]
                csn_start = -1
                csn_end = -1
                if b_rc == True:
                    csn_start = int(fields[13])
                    csn_end = int(fields[12])
                else:
                    csn_start = int(fields[11])
                    csn_end = int(fields[12])

                if abs(end_pos-start_pos)<i_min_len:
                    continue

                if chrm not in self.m_rmsk_annotation:
                    self.m_rmsk_annotation[chrm] = {}
                if start_pos in self.m_rmsk_annotation[chrm]:  ##this is not allowed!!!!
                    print(("Position {0}:{1} has more than  1 annotation!!!!".format(chrm, start_pos)))

                extd_start_pos = start_pos - i_extnd
                extd_end_pos=end_pos+i_extnd
                if extd_start_pos not in self.m_rmsk_annotation[chrm]:
                    self.m_rmsk_annotation[chrm][extd_start_pos] = []
                self.m_rmsk_annotation[chrm][extd_start_pos].append(
                    (extd_end_pos, b_rc, sub_type, csn_start, csn_end))


    #Index by interval tree
    def index_rmsk_annotation_interval_tree(self):
        for chrm in self.m_rmsk_annotation:
            interval_tree = IntervalTree()
            for pos in self.m_rmsk_annotation[chrm]:
                end_pos=self.m_rmsk_annotation[chrm][pos][0][0]
                interval_tree.addi(pos, end_pos)
            self.m_interval_tree[chrm] = interval_tree

    def is_within_repeat_region_interval_tree(self, chrm1, pos):
        chrm=self._process_chrm_name(chrm1)
        if chrm not in self.m_interval_tree:
            return False, -1
        tmp_tree = self.m_interval_tree[chrm]
        set_rslt = tmp_tree[pos]
        if len(set_rslt)==0:
            return False, -1
        for rcd in set_rslt:
            start_pos=rcd[0]
            return True, start_pos
        return False, -1

