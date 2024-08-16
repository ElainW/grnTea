##09/05/2018
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong.simon.chu@gmail.com

# Two classes are defined in this script:
# XReference and XChromosome

import os
import sys
import pysam
from subprocess import *
from multiprocessing import Pool
import global_values
from cmd_runner import *

####
class XChromosome():
    def is_decoy_contig_chrms(self, chrm):
        fields = chrm.split("_")
        if len(fields) > 1:
            return True
        elif chrm == "hs37d5":
            return True

        if len(chrm)>3 and ((chrm[:3] == "HLA") or (chrm[:3] == "HPV") or (chrm[:3] == "HIV") or (chrm[:3] == "CMV")
                            or (chrm[:3] == "CMV") or (chrm[:3] == "MCV") or (chrm[:2] == "SV") or (chrm[:4] == "KSHV")
                            or (chrm[:5] == "decoy") or (chrm[:6] == "random")):#this is some special fields in the bam file
            return True
        fields=chrm.split("-")
        if len(fields)>1:
            return True

        if ("hs38" in chrm) or ("hs37" in chrm):
            return True
####
        # if this is not to call mitochondrial insertion, then filter out chrm related reads
        if global_values.GLOBAL_MITCHONDRION_SWITCH=="OFF":
            if chrm=="MT" or chrm=="chrMT" or chrm=="chrM":#doesn't consider the mitchrondrial DNA
                #print "[TEST]: global value is off"
                return True

        dot_fields = chrm.split(".")
        if len(dot_fields) > 1:
            return True
        else:
            return False


class XReference():
    def __init__(self):
        self.cmd_runner = CMD_RUNNER()

    ## "self.b_with_chr" is the format gotten from the alignment file
    ## all other format should be changed to consistent with the "b_with_chr"
    def process_chrm_name(self, chrm, b_with_chr):
        b_chrm_with_chr = False
        if len(chrm) > 3 and chrm[:3] == "chr":  ##Here remove the "chr"
            b_chrm_with_chr = True

        # print chrm, self.b_with_chr, b_chrm_with_chr ###############################################################
        if b_with_chr == True and b_chrm_with_chr == True:
            return chrm
        elif b_with_chr == True and b_chrm_with_chr == False:
            return "chr" + chrm
        elif b_with_chr == False and b_chrm_with_chr == True:
            return chrm[3:]
        else:
            return chrm

####
    #given chromosome lengths, break to bins
    def break_ref_to_bins(self, m_chrm_length, bin_size):
        m_bins={}
        for chrm in m_chrm_length:
            if chrm not in m_bins:
                m_bins[chrm]=[]
            chrm_lenth=m_chrm_length[chrm]
            tmp_pos=0
            while tmp_pos<chrm_lenth:
                block_end=tmp_pos+bin_size
                if block_end>chrm_lenth:
                    block_end=chrm_lenth
                m_bins[chrm].append((tmp_pos, block_end))
                tmp_pos=block_end
        return m_bins

    ##get index by chrm position
    def get_bin_by_pos(self, chrm, pos, bin_size, m_bins):
        bin_index=-1
        if chrm not in m_bins:
            return bin_index
        bin_index = pos/bin_size
        return m_bins[chrm][bin_index] #return [start, end) of the bin
####
