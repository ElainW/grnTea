#!/bin/python3
# 2021/05/10
##@@author: Yilan (Elain) Wang, Harvard University
##@@contact: yilanwang@g.harvard.edu
'''
1. modified from 2.0count_clip_disc_polyA.py, put everything into a python class
2. add counts of low MAPQ ratio (default MAPQ <= 5)
'''
# import python standard libraries
import os
import sys
import pysam
from multiprocessing import Pool
from subprocess import *
import numpy as np
import re
from shutil import rmtree

# import modules
import global_values
from cmd_runner import *


def unwrap_extract_features_by_chr(arg, **kwarg):
    return Feature_Matrix.extract_features_by_chr(*arg, **kwarg)


class Feature_Matrix():
    
    def __init__(self, input, output, wfolder, bam_list, ref, n_jobs=8, b_train=False, margin=50, cMAPQ=12, low_MAPQ=5, err_margin=15, check_polyA_seq_max=20): # YW 2021/11/2 added b_train
        self.input = input
        self.output = output
        self.wfolder = wfolder
        self.bam_list = bam_list
        self.ref = ref
        self.n_jobs = n_jobs
        self.b_train = b_train
        self.margin = margin
        self.cMAPQ = cMAPQ
        self.low_MAPQ = low_MAPQ
        self.err_margin = err_margin
        self.check_polyA_seq_max = check_polyA_seq_max
        self.cmd_runner = CMD_RUNNER() # to run command from the command line
        self.disc_dict = {} # stores information from after disc info collection and clip-cns + disc-cns filtering
    
    
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
    
    
    # YW 2021/09/24 adapted from x_intermediate_sites.py function call_peak_candidate_sites_calc_std_deviation
    def calculate_clip_pos_std(self, chrm_original_sites, start_pos, end_pos, chrm):
        l_pos = []
        for pos in range(start_pos, end_pos):
            if pos in chrm_original_sites:
                l_pos.extend([pos] * (chrm_original_sites[pos][0] + chrm_original_sites[pos][1]))
        b = np.array(l_pos)
        # print(f"{chrm}:{start_pos}-{end_pos}")
        # print(b)
        f_std = round(np.std(b), 2)
        return f_std
    
        
    def extract_features_by_chr(self, record):
        bam = record[0]
        chrm = record[1]
        swfolder = record[2]
        chrm_original_sites = record[3]
        
        extra_features = {}
    
        samfile = pysam.AlignmentFile(bam, "rb", reference_filename=self.ref)
        for insertion_pos in self.disc_dict[chrm]:
            start_pos = insertion_pos - self.margin
            if start_pos <= 0:
                start_pos = 1
            end_pos = insertion_pos + self.margin + 1
            # initialize the output counts and lists
            soft_clip_len = []
            n_low_MAPQ = 0
            n_total = 0
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
                n_total += 1
                if algnmt.mapping_quality <= self.low_MAPQ:
                    n_low_MAPQ += 1
                    continue
                elif algnmt.mapping_quality < self.cMAPQ:
                    continue
                
                # avoid double count polyA
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
                    c_start_pos = -1 * l_cigar[-1][1]
                    clipped_seq = query_seq[c_start_pos:]
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
            
            # calculate clip_pos_std +/- margin with chrm_original_sites
            clip_pos_std = self.calculate_clip_pos_std(chrm_original_sites, start_pos, end_pos, chrm)
            
            # add the output to the final dictionary
            if len(self.disc_dict[chrm][insertion_pos]) == 21: # YW 2021/09/24 change from 22
                extra_features[insertion_pos] = list(map(str, [longest_soft_clip_len, l_cov, r_cov, polyA, l_polyA, r_polyA, clip_pos_std, n_low_MAPQ, n_total])) # YW 2021/09/24 added clip_pos_std
                # self.disc_dict[chrm][insertion_pos].extend([int(longest_soft_clip_len), l_cov, r_cov, polyA, l_polyA, r_polyA])
            else:
                print(self.disc_dict[chrm][insertion_pos])
                sys.exit(f"The number of fields isn't correct. It is {len(self.disc_dict[chrm][insertion_pos])}!")
                
        samfile.close()
        
        chrm_feat = swfolder + chrm + "_feature_matrix.txt"
        with open(chrm_feat, 'w') as fout:
            for insertion_pos in extra_features:
                fout.write("\t".join([chrm, str(insertion_pos), ""]))
                fout.write("\t".join(extra_features[insertion_pos]) + "\n")
        extra_features.clear() # release memory
    
    
    def extract_features_of_given_list(self, bam, swfolder, m_original_sites):
        '''
        operate on each bam file
        2021/11/03 this function cannot handle very long lists of l_chrm_records
        '''
        l_chrm_records = []
    
        for chrm in self.disc_dict:
            l_chrm_records.append((bam, chrm, swfolder, m_original_sites[chrm]))
    
        pool = Pool(self.n_jobs)
        pool.map(unwrap_extract_features_by_chr, list(zip([self] * len(l_chrm_records), l_chrm_records)), 1)
        pool.close()
        pool.join()
        
        with open(swfolder + "feature_matrix.txt", 'w') as fout:
            for chrm in self.disc_dict:
                chrm_feat = swfolder + chrm + "_feature_matrix.txt"
                with open(chrm_feat, 'r') as fin:
                    for line in fin:
                        fout.write(line)
                os.remove(chrm_feat)
    
    
    def load_in_candidate_list(self):
        '''
        2021/11/03 this function can lead to large memory consumption by storing all chr in one dict
        solution 1: can split the input file into several equal chunks (by line count); make sure to change open('w') to open('a')!!!
        solution 2: can split each chromosome into several equal chunks (by bin_size); make sure to change open('w') to open('a')!!!
        '''
        with open(self.input) as fin_candidate_sites:
            for line in fin_candidate_sites:
                if line[0] == "#": # skip the header line
                    continue
                fields = line.split()
                if len(fields) < 24 and self.b_train: # YW 2021/09/24 change from 25
                    print(fields, " does not have enough fields")
                    continue
                if len(fields) < 23 and not self.b_train: # YW 2021/11/2 added this if statement
                    print(fields, " does not have enough fields")
                    continue
                chrm = fields[0]
                pos = int(fields[1])
                if chrm not in self.disc_dict:
                    self.disc_dict[chrm] = {}
                if pos not in self.disc_dict[chrm]:
                    if self.b_train:
                        self.disc_dict[chrm][pos] = fields[3:]
                    else:
                        self.disc_dict[chrm][pos] = fields[2:]
                
    
    
    def output_sample_feature_matrix_multibam(self, cnt):
        '''
        operate on all bam files in the bam list, i.e., there can be more than 1 bam files
        '''
        # write the feature matrix header
        colnames = "\t".join(["#chr", "pos", "pos+1", "lclip", "rclip", "cr_Alu", "cr_L1", "cr_SVA", "cns_Alu", "cns_L1", "cns_SVA", "raw_ldisc", "raw_rdisc", "ldisc_Alu", "rdisc_Alu", "ldisc_L1", "rdisc_L1", "ldisc_SVA", "rdisc_SVA", "ratio_lcluster", "ratio_rcluster", "dr_Alu", "dr_L1", "dr_SVA", "longest_clip_len", "l_cov", "r_cov", "polyA", "cns_std_l_polyA", "cns_std_r_polyA", "clip_pos_std", "ratio_low_MAPQ"])
        with open(self.output, 'w') as fout:
            fout.write(colnames + "\n")
            # load extra features for each sample into self.disc_dict
            for i in range(cnt):
                swfolder = self.wfolder + str(i) + "/"
                sample_feat_mat = swfolder + "feature_matrix.txt"
                with open(sample_feat_mat, 'r') as fin:
                    for line in fin:
                        fields = line.rstrip().split("\t")
                        chrm = fields[0]
                        insertion_pos = int(fields[1])
                        longest_soft_clip_len = int(fields[2])
                        l_cov = float(fields[3])
                        r_cov = float(fields[4])
                        polyA = int(fields[5])
                        l_polyA = int(fields[6])
                        r_polyA = int(fields[7])
                        clip_pos_std = float(fields[8])
                        n_low_MAPQ = int(fields[9])
                        n_total = int(fields[10])
                        
                        if chrm in self.disc_dict:
                            if insertion_pos in self.disc_dict[chrm]:
                                if len(self.disc_dict[chrm][insertion_pos]) == 21: # change from 22
                                    self.disc_dict[chrm][insertion_pos].extend([longest_soft_clip_len, l_cov, r_cov, polyA, l_polyA, r_polyA, clip_pos_std, n_low_MAPQ, n_total])
                                elif len(self.disc_dict[chrm][insertion_pos]) == 29: # change from 30
                                    self.disc_dict[chrm][insertion_pos][-9] = max(longest_soft_clip_len, self.disc_dict[chrm][insertion_pos][-8])
                                    self.disc_dict[chrm][insertion_pos][-8] += l_cov
                                    self.disc_dict[chrm][insertion_pos][-7] += r_cov
                                    self.disc_dict[chrm][insertion_pos][-6] += polyA
                                    self.disc_dict[chrm][insertion_pos][-5] += l_polyA
                                    self.disc_dict[chrm][insertion_pos][-4] += r_polyA
                                    self.disc_dict[chrm][insertion_pos][-3] = max(clip_pos_std, self.disc_dict[chrm][insertion_pos][-3])
                                    self.disc_dict[chrm][insertion_pos][-2] += n_low_MAPQ
                                    self.disc_dict[chrm][insertion_pos][-1] += n_total
                            else:
                                sys.exit(f"{insertion_pos} does not exist in {chrm} of the disc input")
                        else:
                            sys.exit(f"{chrm} does not exist in the disc input")
                rmtree(swfolder)
            
            # write self.disc_dict into final output
            for chrm in self.disc_dict:
                for insertion_pos in self.disc_dict[chrm]:
                    fout.write("\t".join([chrm, str(insertion_pos), str(insertion_pos + 1), ""]))
                    try:
                        low_MAPQ_ratio = str(self.disc_dict[chrm][insertion_pos][-2]/self.disc_dict[chrm][insertion_pos][-1])
                    except ZeroDivisionError:
                        low_MAPQ_ratio = '-1'
                    fout.write("\t".join(map(str, self.disc_dict[chrm][insertion_pos][:-2])) + "\t" + low_MAPQ_ratio + "\n")
    
    
    def output_sample_feature_matrix(self, cnt=1):
        '''
        operate on all bam files in the bam list, assuming cnt=1
        '''
        # write the feature matrix header
        colnames = "\t".join(["#chr", "pos", "pos+1", "lclip", "rclip", "cr_Alu", "cr_L1", "cr_SVA", "cns_Alu", "cns_L1", "cns_SVA", "raw_ldisc", "raw_rdisc", "ldisc_Alu", "rdisc_Alu", "ldisc_L1", "rdisc_L1", "ldisc_SVA", "rdisc_SVA", "ratio_lcluster", "ratio_rcluster", "dr_Alu", "dr_L1", "dr_SVA", "longest_clip_len", "l_cov", "r_cov", "polyA", "cns_std_l_polyA", "cns_std_r_polyA", "clip_pos_std", "ratio_low_MAPQ"])
        swfolder = self.wfolder + "0/"
        sample_feat_mat = swfolder + "feature_matrix.txt"
        with open(self.output, 'w') as fout, open(sample_feat_mat, 'r') as fin:
            fout.write(colnames + "\n")
            for line in fin:
                fields = line.rstrip().split("\t")
                chrm = fields[0]
                insertion_pos = int(fields[1])
                longest_soft_clip_len = int(fields[2])
                l_cov = float(fields[3])
                r_cov = float(fields[4])
                polyA = int(fields[5])
                l_polyA = int(fields[6])
                r_polyA = int(fields[7])
                clip_pos_std = float(fields[8])
                n_low_MAPQ = int(fields[9])
                n_total = int(fields[10])
                try:
                    low_MAPQ_ratio = str(n_low_MAPQ/n_total)
                except ZeroDivisionError:
                    low_MAPQ_ratio = '-1'
                
                if chrm in self.disc_dict:
                    if insertion_pos in self.disc_dict[chrm]:
                        fout.write("\t".join([chrm, str(insertion_pos), str(insertion_pos + 1), ""]))
                        fout.write("\t".join(self.disc_dict[chrm][insertion_pos]) + "\t")
                        fout.write("\t".join([longest_soft_clip_len, l_cov, r_cov, polyA, l_polyA, r_polyA, clip_pos_std, low_MAPQ_ratio]) + "\n")
                    else:
                        sys.exit(f"{insertion_pos} does not exist in {chrm} of the disc input")
                else:
                    sys.exit(f"{chrm} does not exist in the disc input")
        rmtree(swfolder)

    
    # YW 2021/09/24 added m_original_sites to calculate clip_pos_std (different definition as in x_intermediate_sites.py) here
    def run_feature_extraction(self, m_original_sites):
        print(f"With minimal MAPQ required for clipped reads inspection: {self.cMAPQ}; MAPQ < {self.low_MAPQ} for low MAPQ reads; maximum distance {self.err_margin} bp between clipping position and breakpoint allowed")
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
                self.extract_features_of_given_list(bam, swfolder, m_original_sites)
                cnt += 1
        if cnt == 1:
            self.output_sample_feature_matrix()
        else:
            self.output_sample_feature_matrix_multibam(cnt)