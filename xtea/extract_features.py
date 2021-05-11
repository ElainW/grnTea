#!/bin/python3
# 2021/05/10
# modified from 2.0count_clip_disc_polyA.py, put everything into a python class
# import python standard libraries
import os
import sys
import pysam
import argparse
from multiprocessing import Pool
from subprocess import *
import numpy as np
import re
# import modules
import global_values
from cmd_runner import *


def unwrap_extract_features_by_chr(arg, **kwarg):
    return Feature_Matrix.extract_features_by_chr(*arg, **kwarg)


class Feature_Matrix():
    
    def __init__(self, input, output, wfolder, bam_list, ref, n_jobs=8, margin=50, cMAPQ=12, err_margin=15, check_polyA_seq_max=20):
        self.input = input
        self.output = output
        self.wfolder = wfolder
        self.bam_list = bam_list # integrate this with xTea scripts
        self.ref = ref
        self.n_jobs = n_jobs
        self.margin = margin
        self.cMAPQ = cMAPQ
        self.err_margin = err_margin
        self.check_polyA_seq_max = check_polyA_seq_max
        self.cmd_runner = CMD_RUNNER() # to run command from the command line
        self.disc_dict = {}
    
    
    def is_poly_A_T(self, seq):  ###for a temp version here
        cnt_A = 0
        cnt_T = 0
        mask = 0
        if len(seq) < 4:
            return False
        for ch in seq:
            if ch == 'A' or ch == 'a':
                cnt_A += 1
            elif ch == 'T' or ch == 't':
                cnt_T += 1
            elif ch == 'N':
                mask += 1
        cnt_cutoff = (len(seq)-mask) * 0.75
        if cnt_A + cnt_T > cnt_cutoff and len(seq)>=4:
            return True
        return False
    
    
    # added 2020/08/20 to see how many gold std set have polyA by cns filtering standards
    def is_consecutive_polyA_new(self, seq, b_clip_part_rc, len_seq): 
        if len_seq >= 5:
            # if ("AAAAA" in seq) or ("TAAAA" in seq) or ("ATAAA" in seq) or ("AATAA" in seq) or ("AAATA" in seq) or ("AAAAT" in seq):
            if re.search(r"[ATCG]AAAA", seq) or re.search(r"A[TCG]AAA", seq) or re.search(r"AA[TCG]AA", seq) or re.search(r"AAA[TCG]A", seq) or re.search(r"AAAA[TCG]", seq) or re.search(r"[ATCG]TTTT", seq) or re.search(r"T[ACG]TTT", seq) or re.search(r"TT[ACG]TT", seq) or re.search(r"TTT[ACG]T", seq) or re.search(r"TTTT[ACG]", seq):
                return True
        elif len_seq >= 3:
            if re.search(rf"A{len_seq}", seq) or re.search(rf"T{len_seq}", seq):
                return True
        return False
    
    
    def calc_clip_pos(self, l_cigar, lrclip, map_start):
        clip_pos = 0
        offset = 0
        if lrclip=="L":
            clip_pos = map_start
        elif lrclip=="R":
            for (type, lenth) in l_cigar[:-1]:
                if type == 4 or type == 5 or type == 1:  # (1 for insertion)
                    continue
                else:
                    offset += lenth
            clip_pos = map_start + offset
        return clip_pos
    
    
    def _get_map_end_pos(self, l_cigar, start_pos):
        '''
        function copied from x_coverage.py
        '''
        end_pos=start_pos
        for (type, lenth) in l_cigar:
            if type == 4 or type == 5 or type == 1:  # (1 for insertion)
                continue
            else:
                end_pos += lenth
        return end_pos


    def compute_coverage(self, mcov, win_size):
        '''
        function copied from x_coverage.py
        '''
        isum=0
        for tmp_pos in mcov:
            isum+=mcov[tmp_pos]
        fcov=float(isum)/float(win_size)
        return fcov
    
    
    def cal_coverage(self, l_cigar, map_pos, lcov, rcov, l_start, l_end, r_start, r_end):
        '''
        function modified from x_coverage.py _calc_depth_one_site
        '''
        tmp_end_pos=self._get_map_end_pos(l_cigar, map_pos)
        for i in range(map_pos, tmp_end_pos):
            if i>=l_start and i<=l_end:
                if i not in lcov:
                    lcov[i]=1
                else:
                    lcov[i]+=1
            if i>=r_start and i<=r_end:
                if i not in rcov:
                    rcov[i]=1
                else:
                    rcov[i]+=1
    
    
    def is_clip_qualified(self, clip_pos, insertion_pos):
        if abs(clip_pos-insertion_pos) <= self.err_margin:
            return True
        else:
            return False
    
        
    def extract_features_by_chr(self, record):
        bam = record[0]
        chrm = record[1]
        swfolder = record[2]
    
        samfile = pysam.AlignmentFile(bam, "rb", reference_filename=self.ref)
        for insertion_pos in self.disc_dict[chrm]:
            start_pos = insertion_pos - self.margin
            if start_pos <= 0:
                start_pos = 1
            end_pos = insertion_pos + self.margin + 1
            # initialize the output counts and lists
            soft_clip_len = []
            l_cov = 0.0
            r_cov = 0.0
            polyA = 0
            l_polyA = 0
            r_polyA = 0
            lrclip = ""
            lcov = {}
            rcov = {}
            for algnmt in samfile.fetch(chrm, start_pos, end_pos):  ##fetch reads mapped to "chrm:start_pos-end_pos"
                len_clip_seq = 0
                clipped_seq = ""
                clip_pos = 0
                if algnmt.is_secondary or algnmt.is_supplementary or algnmt.is_unmapped or algnmt.is_duplicate:
                    continue
    
                l_cigar = algnmt.cigar
                if len(l_cigar) < 1:  # wrong alignment
                    continue
                if algnmt.mapping_quality < self.cMAPQ:
                    continue
                
                # avoid double polyA
                polyA_counted = False
    
                query_seq = algnmt.query_sequence
                query_quality = algnmt.query_qualities
                map_pos = algnmt.reference_start
                is_rc = False
                if algnmt.is_reverse == True:  # is reverse complementary
                    is_rc = True
                
                # 2021/05/10 count left/right local coverage
                self.cal_coverage(l_cigar, map_pos, lcov, rcov, start_pos, insertion_pos, insertion_pos + 1, end_pos)
    
                if l_cigar[0][0] == 4:  #left soft clipped
                    lrclip = "L"
                    clip_pos = self.calc_clip_pos(l_cigar, lrclip, map_pos)
                    if self.is_clip_qualified(clip_pos, insertion_pos) == False:
                        continue
                    clipped_seq = query_seq[:l_cigar[0][1]]
                    if self.is_poly_A_T(clipped_seq):
                        polyA += 1
                        polyA_counted = True
                    len_clip_seq=len(clipped_seq)
                    # 2020/08/09 added
                    for ch in clipped_seq:
                        if ch=="N":
                            len_clip_seq -= 1
                    soft_clip_len.append(len_clip_seq)
    
                if l_cigar[-1][0] == 4:  # right soft clipped
                    ##calculate the exact clip position
                    lrclip = "R"
                    clip_pos = self.calc_clip_pos(l_cigar, lrclip, map_pos)
                    if self.is_clip_qualified(clip_pos, insertion_pos) == False:
                        continue
                    start_pos = -1 * l_cigar[-1][1]
                    clipped_seq = query_seq[start_pos:]
                    if self.is_poly_A_T(clipped_seq) and polyA_counted == False:
                        polyA += 1
                    len_clip_seq = len(clipped_seq)
                    # 2020/08/09 added
                    for ch in clipped_seq:
                        if ch=="N":
                            len_clip_seq -= 1
                    soft_clip_len.append(len_clip_seq)
                
                # 2020/08/20 added to look at polyA by cns filtering standards
                # output is problematic... should check again
                s_clip_seq_ck = ""
                if is_rc==False:#not reverse complementary
                    if lrclip == "R":  # right clip
                        if len_clip_seq > self.check_polyA_seq_max:
                            s_clip_seq_ck = clipped_seq[:self.check_polyA_seq_max]
                        elif len_clip_seq > 0:
                            s_clip_seq_ck = clipped_seq    
                    elif lrclip == "L":  # left clip
                        if len_clip_seq > self.check_polyA_seq_max:
                            s_clip_seq_ck = clipped_seq[-1 * self.check_polyA_seq_max:]
                        elif len_clip_seq > 0:
                            s_clip_seq_ck = clipped_seq
                else:# clip part is reverse complementary
                    if lrclip == "L":  # left clip
                        if len_clip_seq > self.check_polyA_seq_max:
                            s_clip_seq_ck = clipped_seq[:self.check_polyA_seq_max]
                        elif len_clip_seq > 0:
                            s_clip_seq_ck = clipped_seq
                    elif lrclip == "R":  # right clip
                        if len_clip_seq > self.check_polyA_seq_max:
                            s_clip_seq_ck = clipped_seq[-1 * self.check_polyA_seq_max:]
                        elif len_clip_seq > 0:
                            s_clip_seq_ck = clipped_seq
                b_polya = self.is_consecutive_polyA_new(s_clip_seq_ck, is_rc, len_clip_seq)
                if b_polya and self.is_clip_qualified(clip_pos, insertion_pos):
                    if lrclip =="L":
                        l_polyA += 1
                    elif lrclip == "R":
                        r_polyA += 1
            
            l_cov = self.compute_coverage(lcov, self.margin)
            r_cov = self.compute_coverage(rcov, self.margin)
            longest_soft_clip_len = 0
            if len(soft_clip_len) > 0:
                longest_soft_clip_len = np.amax(np.array(soft_clip_len))
                # print(f"{longest_soft_clip_len}, type is {type(longest_soft_clip_len)}")
            
            # add the output to the final dictionary
            if len(self.disc_dict[chrm][insertion_pos]) == 22:
                self.disc_dict[chrm][insertion_pos].extend([str(longest_soft_clip_len), str(l_cov), str(r_cov), str(polyA), str(l_polyA), str(r_polyA)])
                # self.disc_dict[chrm][insertion_pos].extend([int(longest_soft_clip_len), l_cov, r_cov, polyA, l_polyA, r_polyA])
            # elif len(self.disc_dict[chrm][insertion_pos]) == 28:
            #     self.disc_dict[chrm][insertion_pos][-6] = max(int(longest_soft_clip_len), disc_chrm_dict[insertion_pos][-6])
            #     self.disc_dict[chrm][insertion_pos][-5] += l_cov
            #     self.disc_dict[chrm][insertion_pos][-4] += r_cov
            #     self.disc_dict[chrm][insertion_pos][-3] += polyA
            #     self.disc_dict[chrm][insertion_pos][-2] += l_polyA
            #     self.disc_dict[chrm][insertion_pos][-1] += r_polyA
            else:
                sys.exit(f"The number of fields isn't correct. It is {len(disc_chrm_dict[insertion_pos])}!")
                
        samfile.close()
        
        chrm_feat = swfolder + chrm + "_feature_matrix.txt"
        with open(chrm_feat, 'w') as fout:
            for insertion_pos in self.disc_dict[chrm]:
                fout.write("\t".join([chrm, str(insertion_pos), ""]))
                fout.write("\t".join(self.disc_dict[chrm][insertion_pos]) + "\n")
    
    
    def extract_features_of_given_list(self, bam, swfolder):
        l_chrm_records = []
    
        for chrm in self.disc_dict:
            l_chrm_records.append((bam, chrm, swfolder))
    
        pool = Pool(self.n_jobs)
        pool.map(unwrap_extract_features_by_chr, list(zip([self] * len(l_chrm_records), l_chrm_records)), 1)
        pool.close()
        pool.join()
    
    
    def load_in_candidate_list(self):
        with open(self.input) as fin_candidate_sites:
            for line in fin_candidate_sites:
                if line[0] == "#": # skip the header line
                    continue
                fields = line.split()
                if len(fields) < 24:
                    print(fields, " does not have enough fields")
                    continue
                chrm = fields[0]
                pos = int(fields[1])
                if chrm not in self.disc_dict:
                    self.disc_dict[chrm] = {}
                if pos not in self.disc_dict[chrm]:
                    self.disc_dict[chrm][pos] = fields[2:]
    
    
    def output_sample_feature_matrix(self, swfolder):
        # write the feature matrix header
        colnames = "\t".join(["#chr", "pos", "lclip", "rclip", "cr_Alu", "cr_L1", "cr_SVA", "cns_Alu", "cns_L1", "cns_SVA", "clip_pos_std", "raw_ldisc", "raw_rdisc", "ldisc_Alu", "rdisc_Alu", "ldisc_L1", "rdisc_L1", "ldisc_SVA", "rdisc_SVA", "ratio_lcluster", "ratio_rcluster", "dr_Alu", "dr_L1", "dr_SVA", "longest_clip_len", "l_cov", "r_cov", "polyA", "cns_std_l_polyA", "cns_std_r_polyA"])
        with open(swfolder + "feature_matrix.txt", 'w') as fout:
            fout.write(colnames + "\n")
            for chrm in self.disc_dict:
                chrm_feat = swfolder + chrm + "_feature_matrix.txt"
                with open(chrm_feat, 'r') as fin:
                    for line in fin:
                        fout.write(line)
                    os.remove(chrm_feat)


    def run_feature_extraction(self):
        '''
        now cannot handle the case where there are multiple bam files for one sample
        '''
        print(f"With minimal MAPQ required for clipped reads inspection: {self.cMAPQ}; maximum distance {self.err_margin} bp between clipping position and breakpoint allowed")
        self.load_in_candidate_list()
        cnt = 0
        with open(self.bam_list, 'r') as b_list:
            for lines in b_list:
                bam, read_type = lines.rstrip().split('\t')
                if read_type != "illumina":
                    raise NotImplementedError
                swfolder = self.wfolder + str(cnt) + "/"
                scmd = f"mkdir -p {swfolder}"
                self.cmd_runner.run_cmd_small_output(scmd)
                self.extract_features_of_given_list(bam, swfolder)
                self.output_sample_feature_matrix(swfolder)