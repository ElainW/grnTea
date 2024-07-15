##11/22/2017
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong.simon.chu@gmail.com

# YW note, classes defined in this script: BamInfo
# 2020/08/02 github updated cnt_discordant_pairs function

import os
import sys
import pysam
from x_intermediate_sites import * # YW 2020/07/03 note: didn't see XIntermediateSites class used here 
from x_reference import *
import global_values
from disc_cluster import *


def unwrap_self_extract_mate_reads_for_region_short(arg, **kwarg):
    return BamInfo.extract_mate_reads_of_region_short(*arg, **kwarg)

####Function: Functions for processing alignment, like trim, change format, and etc.
class BamInfo():#
    def __init__(self, sf_bam, sf_ref):
        self.sf_bam = sf_bam
        self.out_header = None
        self.chrm_id_name = {}
        self.sf_reference=sf_ref
        # if the chromosome name in format: chr1, then return true,

    def get_all_reference_names(self):
        bamfile = pysam.AlignmentFile(self.sf_bam, "rb", reference_filename=self.sf_reference)
        header = bamfile.header
        bamfile.close()

        l_chrms = header['SQ']
        m_chrms = {}
        for record in l_chrms:
            chrm_name = record['SN']
            m_chrms[chrm_name] = 1
        return m_chrms

    # get the chr name and length in a dictionary
    def get_all_chrom_name_length(self):
        bamfile = pysam.AlignmentFile(self.sf_bam, "rb", reference_filename=self.sf_reference)
        header = bamfile.header
        bamfile.close()

        l_chrms = header['SQ']
        m_chrms = {}
        for record in l_chrms:
            chrm_name = record['SN']
            chrm_length = int(record['LN'])
            m_chrms[chrm_name] = chrm_length
        return m_chrms

####
    # else if in format: 1, then return false
    def is_chrm_contain_chr(self):
        samfile = pysam.AlignmentFile(self.sf_bam, "rb", reference_filename=self.sf_reference)
        self.out_header = samfile.header  # it's a dictionary
        samfile.close()

        l_chrms = self.out_header['SQ']
        m_chrms = {}
        for record in l_chrms:
            chrm_id = record['SN']
            m_chrms[chrm_id] = 1

        b_with_chr = False
        if "chr1" in m_chrms:
            b_with_chr = True

        self.b_with_chr = b_with_chr
        return b_with_chr
####
    ###Note: hard code here, use mapping quality 20 as cutoff for anchor reads
    # YW 2020/08/03 github update
    # YW 2021/04/08 update xannotation input to count the number of mates mapping to repeatmasker annotation coordinates
    # YW 2021/04/16 change output of the function
    # YW 2021/04/26 add collect discordant read pos
    # YW 2021/04/29 add L_R input to indicate which side the of the insertion position we're dealing with
    def cnt_collect_discordant_pairs(self, bamfile, m_chrm_ids, chrm, start, end, i_is, f_dev, xannotation_Alu, xannotation_L1, xannotation_SVA, sf_disc_pos, L_R):
        n_cnt_Alu, n_cnt_L1, n_cnt_SVA = 0, 0, 0
        n_raw_cnt=0 #this save the raw discordant PE pairs, including: abnormal insert size, abnormal direction
        iter_alignmts = bamfile.fetch(chrm, start, end)
        xchrom = XChromosome()
        i_max_is=500 # YW 2020/08/02 changed from 2500
        if int(i_is+3*f_dev)>i_max_is:
            i_max_is=int(i_is+3*f_dev)

        m_mate_pos={}
        with open(sf_disc_pos, 'a') as f_disc_pos: # YW 2021/05/10 changed from 'w' to 'a'
            for algnmt in iter_alignmts:
                if algnmt.is_duplicate == True or algnmt.is_supplementary == True:  ##skip duplicate and supplementary ones
                    continue
                if algnmt.is_unmapped == True or algnmt.mate_is_unmapped == True:  #### for now, just skip the unmapped reads
                    continue
                if algnmt.is_secondary == True:  ##skip secondary alignment
                    continue
                if algnmt.mapping_quality < global_values.MINIMUM_DISC_MAPQ:###############anchor mapping quality
                    continue
                if algnmt.next_reference_id<0:
                    continue
    ####
                map_pos=algnmt.reference_start
                if algnmt.next_reference_id not in m_chrm_ids:
                    continue
                mate_chrm = algnmt.next_reference_name
                if xchrom.is_decoy_contig_chrms(mate_chrm) == True:
                    continue
                mate_pos = algnmt.next_reference_start
    
    
                #####This version only count the number of discordant pairs whose two reads aligned to different chroms
                #####Or on the same chromosome, but aligned quite far away (by default >1M) => YW 2021/04/30 > 500 bp by default
                #Note, here we use "start" to represent the "map_pos" of each read.
                if (chrm != mate_chrm) or (chrm == mate_chrm and abs(mate_pos - start) > max(i_is, i_max_is)):
                    if mate_chrm not in m_mate_pos:
                        m_mate_pos[mate_chrm]=[]
                    m_mate_pos[mate_chrm].append(mate_pos)
                    n_raw_cnt+=1
                    b_mate_within_rep_Alu, rep_start_mate_Alu = xannotation_Alu.is_within_repeat_region_interval_tree(mate_chrm, mate_pos)
                    b_mate_within_rep_L1, rep_start_mate_L1 = xannotation_L1.is_within_repeat_region_interval_tree(mate_chrm, mate_pos)
                    b_mate_within_rep_SVA, rep_start_mate_SVA = xannotation_SVA.is_within_repeat_region_interval_tree(mate_chrm, mate_pos)
                    
                    ##########################################
                    # YW 2021/04/26 added to write to disc_pos
                    query_name = algnmt.query_name
                    b_first = False if algnmt.is_read2 else True
                    s_mate_first = 0 if b_first else 1 # whether the mate read is the "first read" in a pair
                    # YW 2021/04/29 added to distinguish left disc from right disc
                    insertion_pos = end if L_R == "left" else start - 1

                    # here mate_chrm must be the style in the bam file
                    # And chrm must be the style in the candidate file
                    # YW 2021/04/29 take out irrelevant info
                    s_mate_pos_info = "\t".join(list(map(str, [chrm, map_pos, mate_chrm, mate_pos, s_mate_first, query_name, insertion_pos])))
                    f_disc_pos.write(s_mate_pos_info + "\n")
                    ###########################################
                
                    # seq = alignmt.query_sequence
                    if b_mate_within_rep_Alu:
                        n_cnt_Alu += 1
                        # YW TO DO: write out the seq
                    if b_mate_within_rep_L1:
                        n_cnt_L1 += 1
                        # YW TO DO: write out the seq
                    if b_mate_within_rep_SVA:
                        n_cnt_SVA += 1
                        # YW TO DO: write out the seq
        dc = DiscCluster()
        # YW 2021/04/16 changed the following few lines
        # b_cluster, c_chrm, c_pos=dc.form_one_side_cluster(m_mate_pos, i_is, global_values.MIN_RAW_DISC_CLUSTER_RATIO)
        cluster_ratio, c_chrm, c_pos = dc.form_one_side_cluster(m_mate_pos, i_is)
        # return n_cnt, n_raw_cnt, (b_cluster, c_chrm, c_pos)
        return n_raw_cnt, n_cnt_Alu, n_cnt_L1, n_cnt_SVA, (cluster_ratio, c_chrm, c_pos)

    ## "self.b_with_chr" is the format gotten from the alignment file
    ## all other format should be changed to consistent with the "self.b_with_chr"
    def process_chrm_name(self, chrm, b_with_chr):
        b_chrm_with_chr = False
        if len(chrm) > 3 and chrm[:3] == "chr":  ##Here remove the "chr"
            b_chrm_with_chr = True

        # print chrm, self.b_with_chr, b_chrm_with_chr #######################################################################
        if b_with_chr == True and b_chrm_with_chr == True:
            return chrm
        elif b_with_chr == True and b_chrm_with_chr == False:
            return "chr" + chrm
        elif b_with_chr == False and b_chrm_with_chr == True:
            return chrm[3:]
        else:
            return chrm

    
    # YW 2021/04/29 abbreviate this function from _parse_read_name_pos_of_region
    def _parse_read_name_pos_of_region_short(self, sf_name_pos, chrm, start, end):
        m_reads={}
        with open(sf_name_pos) as fin_name_pos:
            # each line in format: chrm, map_pos, mate_chrm, mate_pos, query_name, insertion_pos
            for line in fin_name_pos:
                fields=line.split()
                ins_chrm = fields[0]
                anchor_read_map_pos=int(fields[1]) ##this is the map_pos of the anchor read
                mate_chrm=fields[2]
                mate_pos=int(fields[3])
                
                # YW 2021/04/29 commented out below, because we want to retain disc reads mapping to different chr
                # if chrm != mate_chrm:
                #     continue
                if mate_pos<(start-100) or mate_pos>(end+100):
                    continue
                s_first=fields[4]
                rname=fields[5]
                s_insertion_pos=int(fields[6])

                s_insertion_site = f"{ins_chrm}~{s_insertion_pos}"
                s_info = f"{rname}~{s_first}"
                if s_info not in m_reads:
                    m_reads[s_info]={}
                m_reads[s_info][s_insertion_site]=anchor_read_map_pos
        return m_reads
    
    # YW 2021/04/29 abbreviate the function from extract_mate_reads_of_region
    # YW 2021/04/29 now the read info format is >read_name~is_first~anchor_map_pos~insertion_pos
    def extract_mate_reads_of_region_short(self, record):
        chrm=record[0]
        bin_start=int(record[1])
        bin_end=int(record[2])
        sf_name_pos=record[3]
        s_working_folder=record[4]

        m_reads=self._parse_read_name_pos_of_region_short(sf_name_pos, chrm, bin_start, bin_end)

        #first write the bam to file, and then index it, and finally get the fastq reads
        bamfile = pysam.AlignmentFile(self.sf_bam, "rb", reference_filename=self.sf_reference)
        sf_tmp_disc_fa = s_working_folder + f"{chrm}_{bin_start}_{bin_end}.tmp.disc.fa"
        with open(sf_tmp_disc_fa,"w") as fout_fa:
            for alignmt in bamfile.fetch(chrm, bin_start, bin_end):
                s_read_name = alignmt.query_name
                read_seq = alignmt.query_sequence
                is_first = 1 if alignmt.is_read1 else 0

                s_read_id = f"{s_read_name}~{is_first}"
                if s_read_id not in m_reads:
                    continue
                for s_insertion in m_reads[s_read_id]:
                    anchor_map_pos=m_reads[s_read_id][s_insertion]
                    s_read_info = f">{s_read_id}~{anchor_map_pos}~{s_insertion}\n"
                    fout_fa.write(s_read_info)
                    fout_fa.write(read_seq + "\n")
        bamfile.close()
        
    # YW 2021/04/29 abbreviate the function from extract_mate_reads_by_name
    def extract_mate_reads_by_name_short(self, sf_name_pos, bin_size, sf_working_folder, n_jobs,  sf_out_fa):
        if sf_working_folder[-1]!="/":
            sf_working_folder+="/"
        m_chrm_lenth = self.get_all_chrom_name_length()
        xref = XReference()
        m_bin_pos = xref.break_ref_to_bins(m_chrm_lenth, bin_size) #in format: m{chr:[bin_start, bin_end]}

        xchrom = XChromosome()
        l_bin_info = []
        for chrm in m_bin_pos:
            if xchrom.is_decoy_contig_chrms(chrm)==True:
                continue
            for (bin_start, bin_end) in m_bin_pos[chrm]:
                l_bin_info.append((chrm, bin_start, bin_end, sf_name_pos, sf_working_folder))
        pool = Pool(n_jobs)
        pool.map(unwrap_self_extract_mate_reads_for_region_short, list(zip([self] * len(l_bin_info), l_bin_info)), 1)
        pool.close()
        pool.join()

        #merge the fa files
        with open(sf_out_fa, "w") as fout_disc_fa:
            for record in l_bin_info:
                chrm = record[0]
                istart = record[1]
                iend = record[2]
                sf_tmp_disc_fa = sf_working_folder + f"{chrm}_{istart}_{iend}.tmp.disc.fa"
                with open(sf_tmp_disc_fa) as fin_tmp:
                    for line in fin_tmp:
                        fout_disc_fa.write(line)
                if os.path.isfile(sf_tmp_disc_fa):#clean the temporary file
                    os.remove(sf_tmp_disc_fa)

