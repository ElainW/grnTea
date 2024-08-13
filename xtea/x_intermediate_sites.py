##11/27/2017
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu
'''
# YW notes: defined python class: XIntermediateSites
# YW 2020/08/03 github update:
# added function load_in_candidate_list_one_line, is_in_existing_list function, are_sites_close function
'''
import os
from subprocess import *
import global_values
import numpy as np

# import modules
from cmd_runner import *

class XIntermediateSites():

    # YW 2021/03/19 update how the nearby mate in rep is calculated (now by Alu, L1, SVA)
    def parse_sites_with_clip_cutoff_for_chrm(self, m_clip_pos_freq, cutoff_left_clip, cutoff_right_clip,
                                              cutoff_clip_mate_in_rep, cutoff_clip_mate_in_cns):
        m_candidate_sites = {}
        for pos in m_clip_pos_freq:
            ####here need to check the nearby region
            nearby_left_freq = 0
            nearby_right_freq = 0
            nearby_mate_in_rep_Alu, nearby_mate_in_rep_L1, nearby_mate_in_rep_SVA = 0, 0, 0
            cns_Alu, cns_L1, cns_SVA = 0, 0, 0 # YW 2021/04/23 added
            
            for i in range(-1 * global_values.NEARBY_REGION, global_values.NEARBY_REGION):
                i_tmp_pos = pos + i
                if i_tmp_pos in m_clip_pos_freq:
                    nearby_left_freq += m_clip_pos_freq[i_tmp_pos][0]
                    nearby_right_freq += m_clip_pos_freq[i_tmp_pos][1]
                    nearby_mate_in_rep_Alu += (
                    m_clip_pos_freq[i_tmp_pos][2] + m_clip_pos_freq[i_tmp_pos][5] + m_clip_pos_freq[i_tmp_pos][6])
                    cns_Alu += m_clip_pos_freq[i_tmp_pos][5] + m_clip_pos_freq[i_tmp_pos][6] # YW 2021/04/23 added
                    nearby_mate_in_rep_L1 += (
                    m_clip_pos_freq[i_tmp_pos][3] + m_clip_pos_freq[i_tmp_pos][7] + m_clip_pos_freq[i_tmp_pos][8])
                    cns_L1 += m_clip_pos_freq[i_tmp_pos][7] + m_clip_pos_freq[i_tmp_pos][8] # YW 2021/04/23 added
                    nearby_mate_in_rep_SVA += (
                    m_clip_pos_freq[i_tmp_pos][4] + m_clip_pos_freq[i_tmp_pos][9] + m_clip_pos_freq[i_tmp_pos][10])
                    cns_SVA += m_clip_pos_freq[i_tmp_pos][9] + m_clip_pos_freq[i_tmp_pos][10] # YW 2021/04/23 added

            b_candidate=False
            # if nearby_left_freq >= cutoff_left_clip and nearby_right_freq >= cutoff_right_clip \
            #         and nearby_mate_in_rep >= cutoff_clip_mate_in_rep:
            #     b_candidate=True
            if (nearby_left_freq >= cutoff_left_clip or nearby_right_freq >= cutoff_right_clip) \
                    and (nearby_mate_in_rep_Alu >= cutoff_clip_mate_in_rep or \
                         nearby_mate_in_rep_L1 >= cutoff_clip_mate_in_rep or \
                         nearby_mate_in_rep_SVA >= cutoff_clip_mate_in_rep) \
                    and (cns_Alu >= cutoff_clip_mate_in_cns or \
                         cns_L1 >= cutoff_clip_mate_in_cns or \
                         cns_SVA >= cutoff_clip_mate_in_cns):
                b_candidate=True

            if b_candidate==True:
                # if nearby_left_freq >= cutoff_left_clip and nearby_right_freq >= cutoff_right_clip:
                i_left_cnt = m_clip_pos_freq[pos][0]
                i_right_cnt = m_clip_pos_freq[pos][1]
                i_mate_in_rep_cnt_Alu = m_clip_pos_freq[pos][2]
                i_mate_in_rep_cnt_L1 = m_clip_pos_freq[pos][3]
                i_mate_in_rep_cnt_SVA = m_clip_pos_freq[pos][4]
                m_candidate_sites[pos] = (i_left_cnt, i_right_cnt, i_mate_in_rep_cnt_Alu, i_mate_in_rep_cnt_L1, i_mate_in_rep_cnt_SVA, cns_Alu, cns_L1, cns_SVA)
        return m_candidate_sites
####
    ####

    ###output the candidate list in a file
    def output_candidate_sites(self, m_candidate_list, sf_out):
        with open(sf_out, "w") as fout_candidate_sites:
            for chrm in m_candidate_list:
                if self.is_decoy_contig_chrms(chrm):  ####we are not interested in decoy and other contigs!!!!
                    continue
                for pos in m_candidate_list[chrm]:
                    # lth = len(m_candidate_list[chrm][pos])
                    # fout_candidate_sites.write(chrm + "\t" + str(pos) + "\t")
                    fout_candidate_sites.write("\t".join([chrm, str(pos), ""]))
                    fout_candidate_sites.write("\t".join([str(i) for i in m_candidate_list[chrm][pos]]) + "\n")
                    # for i in range(lth):
                    #     s_feature = str(m_candidate_list[chrm][pos][i])
                    #     fout_candidate_sites.write(s_feature + "\t")
                    # fout_candidate_sites.write("\n")
#
    def is_decoy_contig_chrms(self, chrm):
        fields = chrm.split("_")
        if len(fields) > 1:
            return True
        elif chrm == "hs37d5":
            return True

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


    def load_in_candidate_list(self, sf_candidate_list):
        m_list = {}
        with open(sf_candidate_list) as fin_candidate_sites:
            for line in fin_candidate_sites:
                fields = line.split()
                if len(fields)<10: # 2021/04/01 YW updated from 3, updated from 5
                    print(fields, " does not have enough fields")
                    continue
                chrm = fields[0]
                pos = int(fields[1])
                if chrm not in m_list:
                    m_list[chrm] = {}
                if pos not in m_list[chrm]:
                    # m_list[chrm][pos] = []
                    m_list[chrm][pos] = [int(fields[2]), int(fields[3])] # YW 2021/09/28 changed
                # YW 2021/09/24 changed because only lclip, rclip are used in subsequent steps
                # for ivalue in fields[2:]:
                # for ivalue in fields[2:4]:
                #     m_list[chrm][pos].append(int(ivalue))
        return m_list

    ####
    

    # In the previous step (call_TEI_candidate_sites), some sites close to each other may be introduced together
    # If there are more than 1 site close to each other, than use the peak site as a representative
    def call_peak_candidate_sites(self, m_candidate_sites, peak_window):#
        m_peak_candidate_sites = {}
        for chrm in m_candidate_sites:
            l_pos = list(m_candidate_sites[chrm].keys())
            l_pos.sort()  ###sort the candidate sites
            pre_pos = -1
            set_cluster = set()
            for pos in l_pos:
                if pre_pos == -1:
                    pre_pos = pos
                    set_cluster.add(pre_pos)
                    continue

                if pos - pre_pos > peak_window:  # find the peak in the cluster
                    max_clip = 0
                    tmp_candidate_pos = 0
                    for tmp_pos in set_cluster:
                        tmp_left_clip = int(m_candidate_sites[chrm][tmp_pos][0]) #left clip
                        tmp_right_clip = int(m_candidate_sites[chrm][tmp_pos][1]) #right clip
                        tmp_all_clip = tmp_left_clip + tmp_right_clip #all the clip
                        if max_clip < tmp_all_clip:
                            tmp_candidate_pos = tmp_pos
                            max_clip = tmp_all_clip
                    set_cluster.clear()
                    if chrm not in m_peak_candidate_sites:
                        m_peak_candidate_sites[chrm] = {}
                    if tmp_candidate_pos not in m_peak_candidate_sites[chrm]:
                        m_peak_candidate_sites[chrm][tmp_candidate_pos] = [max_clip]
                pre_pos = pos
                set_cluster.add(pre_pos)

            # push out the last group
            max_clip = 0
            tmp_candidate_pos = 0
            ####check the standard deviation
            for tmp_pos in set_cluster:
                tmp_left_clip = int(m_candidate_sites[chrm][tmp_pos][0])
                tmp_right_clip = int(m_candidate_sites[chrm][tmp_pos][1])
                tmp_all_clip = tmp_left_clip + tmp_right_clip
                if max_clip < tmp_all_clip:
                    tmp_candidate_pos = tmp_pos
                    max_clip = tmp_all_clip
            if chrm not in m_peak_candidate_sites:
                m_peak_candidate_sites[chrm] = {}
            if tmp_candidate_pos not in m_peak_candidate_sites[chrm]:
                ##Here, use list in order to output the list (by call the output_candidate_sites function)
                m_peak_candidate_sites[chrm][tmp_candidate_pos] = [max_clip]
        return m_peak_candidate_sites
####
#
    # YW 2021/04/21 wrote this new function to incorporate info from dict m_sites_clip_peak and sf_raw_disc
    # YW 2021/04/29 added cns_cutoff: clip cns + disc cns >= 1
    # YW 2021/09/24 removed all code related to m_sites_clip_peak
    def merge_clip_disc_new(self, sf_clip, sf_raw_disc, sf_out, sf_out2, cns_cutoff=1): # YW 2021/04/21 alternatively, add info from m_original_sites to m_sites_clip_peak
        # YW 2021/11/16 added to calculate clip_pos_std
        m_clip = dict()
        with open(sf_out, "w") as fout_list:
            # YW 2021/05/10 write the col names
            # YW 2021/09/24 remove "clip_pos_std"
            fout_list.write("\t".join(["#chr", "pos", "lclip", "rclip", "cr_Alu", "cr_L1", "cr_SVA", "cns_Alu", "cns_L1", "cns_SVA", "raw_ldisc", "raw_rdisc", "ldisc_Alu", "rdisc_Alu", "ldisc_L1", "rdisc_L1", "ldisc_SVA", "rdisc_SVA", "ratio_lcluster", "ratio_rcluster", "dr_Alu", "dr_L1", "dr_SVA"]) + "\n")
            m_disc={}
            with open(sf_raw_disc) as fin_disc:
                for line in fin_disc:
                    fields=line.split()
                    # YW 2021/04/21 wrote the following line to accommodate the features in sf_raw_disc:
                    # raw_left_disc, raw_right_disc, s_left_disc_Alu, s_right_disc_Alu, s_left_disc_L1, s_right_disc_L1, s_left_disc_SVA, s_right_disc_SVA, r_lcluster, r_rcluster
                    # YW 2021/04/29 added the following disc reads to cns: dc_Alu, dc_L1, dc_SVA
                    chrm, pos = fields[:2]
                    if chrm not in m_disc:
                        m_disc[chrm]={}
                    m_disc[chrm][pos]=fields[2:]
            with open(sf_clip) as fin_clip:
                # YW 2020/07/20 added this to avoid printing multiple chr doesn't exist message
                chr_not_pre_list = []
                for line in fin_clip:
                    fields = line.split() # YW 2021/04/23: chrm, pos, lclip, rclip, cr_Alu, cr_L1, cr_SVA, cns_Alu, cns_L1, cns_SVA
                    chrm, pos = fields[:2]
                    ###################################################################
                    # YW added 2021/11/16 to calculate clip_pos_std
                    if chrm not in m_clip:
                        m_clip[chrm] = dict()
                    if int(pos) not in m_clip[chrm]:
                        m_clip[chrm][int(pos)] = [int(fields[2]), int(fields[3])]
                    ###################################################################
                    d_fields = ['0', '0', '0', '0', '0', '0', '0', '0', '0.0', '0.0', '0', '0', '0']
                    # YW 2020/07/20 added this to avoid printing multiple times a chr doesn't exist message
                    if chrm not in m_disc:
                        if chrm not in chr_not_pre_list:
                            chr_not_pre_list.append(chrm)
                            print("Error happen at merge clip and disc feature step: {0} not exist. So set left/right disc count to 0".format(chrm))
                        # continue
                    else: 
                        if pos in m_disc[chrm]: # YW 2020/07/20 added and chrm in m_disc
                            d_fields = m_disc[chrm][pos]
                    
                    # YW 2021/04/21 add clip pos std info
                    # YW 2021/09/24 remove clip pos std info here 
                    # f_clip_std = '-1'
                    # if chrm in m_sites_clip_peak:
                    #     if int(pos) in m_sites_clip_peak[chrm]:
                    #         f_clip_std = str(m_sites_clip_peak[chrm][int(pos)][-1])
                    # 
                    # fields.append(f_clip_std)
                    ##########################################
                    # YW 2021/04/29 added this filtering for counts of reads mapping to cns sequence
                    # YW 2021/05/07 now fields[-1] is f_clip_std
                    # YW 2021/09/24 change the numbers because now f_clip_std is removed
                    if int(fields[-3]) + int(d_fields[-3]) >= cns_cutoff or int(fields[-2]) + int(d_fields[-2]) >= cns_cutoff or int(fields[-1]) + int(d_fields[-1]) >= cns_cutoff:
                        fields.extend(d_fields)
                ##########################################
                        fout_list.write("\t".join(fields) + "\n")           
        # YW 2021/09/24 adapted from x_intermediate_sites.py function call_peak_candidate_sites_calc_std_deviation
        def calculate_clip_pos_std(m_clip, f_out): # need to write out a file
            with open(f_out,'w') as fout:
                for chrm in m_clip:
                    for pos in m_clip[chrm]:
                        l_pos = []
                        # print(f"{chrm}:{pos}")
                        for i in range(pos-global_values.MARGIN, pos+global_values.MARGIN+1):
                            if i in m_clip[chrm]:
                                l_pos.extend([i] * (m_clip[chrm][i][0] + m_clip[chrm][i][1]))
                        b = np.array(l_pos)
                        # print(b)
                        f_std = round(np.std(b), 2)
                        fout.write(f"{chrm}\t{pos}\t{f_std}\n") # need to include chrm, pos for correcting sorting and indexing for future df merging
        calculate_clip_pos_std(m_clip, sf_out2)
    
    
    # YW 2021/12/08 added, to add the clip_pos_std info back to the sf_candidate_list
    def add_clip_pos_std(self, f_others, f_clip_pos_std, f_out):
        cmd_runner = CMD_RUNNER()
        sort_f_others = "awk -F'\\t' '{OFS=\"\\t\"; print $1,$2,$2+1,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23}' " + f_others + " | sortBed -i stdin"
        cmd_runner.run_cmd_to_file(sort_f_others, f_others + ".sorted")
        sort_f_clip_pos_std = "awk -F'\\t' '{OFS=\"\\t\"; print $1,$2,$2+1,$3}' " + f_clip_pos_std + " | sortBed -i stdin"
        cmd_runner.run_cmd_to_file(sort_f_clip_pos_std, f_clip_pos_std + ".sorted")
        bedintersect = f"intersectBed -wa -wb -sorted -a {f_others}.sorted -b {f_clip_pos_std}.sorted | cut -f1,2,4-24,28"
        cmd_runner.run_cmd_to_file(bedintersect, f_out)
        os.remove(f_others + ".sorted")
        os.remove(f_clip_pos_std + ".sorted")
        os.remove(f_others)
        os.remove(f_clip_pos_std)

####
