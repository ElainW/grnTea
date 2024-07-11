##11/04/2018
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu
# YW note: defined the following classes: ClipReadInfo, LContigClipReadInfo

import pysam
from x_annotation import *
from x_alignments import *
from x_intermediate_sites import *
import global_values
from x_polyA import *
from bwa_align import *

####
##Function: pool.map doesn't accept function in the class
##So, use this wrapper to solve this
def unwrap_self_cnt_clip_pos(arg, **kwarg):
    return ClipReadInfo.collect_clip_positions_by_chrm(*arg, **kwarg)


# YW 2021/09/29 added with predefined loci
def unwrap_self_collect_clip_info_and_parts(arg, **kwarg):
    return ClipReadInfo.collect_clip_info_and_parts_by_chrm(*arg, **kwarg)


def unwrap_self_collect_clip_parts(arg, **kwarg):
    return ClipReadInfo.collect_clipped_parts_by_chrm(*arg, **kwarg)


def unwrap_self_cnt_clip_parts(arg, **kwarg):
    return ClipReadInfo.run_cnt_clip_part_aligned_to_rep_by_chrm(*arg, **kwarg)


####This class used for collect clip-read related information from the alignment
class ClipReadInfo():
    def __init__(self, sf_bam, n_jobs, sf_ref):
        self.sf_bam = sf_bam
        self.working_folder = "./"
        self.n_jobs = n_jobs
        self.sf_reference = sf_ref

    def set_working_folder(self, working_folder):
        self.working_folder = working_folder
        if working_folder[-1] != "/":
            self.working_folder += "/"


    def _cvt_to_Ascii_quality(self, l_score):
        new_score = [x + 33 for x in l_score]
        return ''.join(map(chr, new_score))

    def _is_qualified_clip(self, l_score):
        n_char=len(l_score)
        n_half=n_char/2
        n_cnt=0
        for i_score in l_score:
            if i_score < global_values.CLIP_PHRED_SCORE_CUTOFF:
                n_cnt+=1
            if n_cnt>n_half:
                return False
        return True


    def _get_chrm_id_name(self, samfile):
        m_chrm = {}
        references = samfile.references
        for schrm in references:
            chrm_id = samfile.get_tid(schrm)
            m_chrm[chrm_id] = schrm
        m_chrm[-1] = "*"
        return m_chrm

    ####return two files:
    ####1. the clip_position
    # YW 2021/03/18 add XAnnotation instances for Alu, L1, SVA
    def collect_clip_positions_by_chrm(self, record):
        chrm = record[0]
        sf_bam = record[1]
        working_folder = record[2]
        sf_annotation_Alu = record[3]
        sf_annotation_L1 = record[4]
        sf_annotation_SVA = record[5]
        b_with_chr = record[6]
        i_clip_cutoff = int(record[7])
        b_se = record[8] #whether this is single end reads
####
        xannotation_Alu = XAnnotation(sf_annotation_Alu)
        xannotation_Alu.set_with_chr(b_with_chr)
        xannotation_L1 = XAnnotation(sf_annotation_L1)
        xannotation_L1.set_with_chr(b_with_chr)
        xannotation_SVA = XAnnotation(sf_annotation_SVA)
        xannotation_SVA.set_with_chr(b_with_chr)
        #xannotation.load_rmsk_annotation()
        i_min_copy_len=0
        i_boundary_extnd=0
        # if global_values.IS_CALL_SVA==True:
        #     i_boundary_extnd=global_values.SVA_ANNOTATION_EXTND
        xannotation_Alu.load_rmsk_annotation_with_extnd_with_lenth_cutoff(i_boundary_extnd, i_min_copy_len)
        xannotation_Alu.index_rmsk_annotation_interval_tree()
        xannotation_L1.load_rmsk_annotation_with_extnd_with_lenth_cutoff(i_boundary_extnd, i_min_copy_len)
        xannotation_L1.index_rmsk_annotation_interval_tree()
        xannotation_SVA.load_rmsk_annotation_with_extnd_with_lenth_cutoff(global_values.SVA_ANNOTATION_EXTND, i_min_copy_len)
        xannotation_SVA.index_rmsk_annotation_interval_tree()
        

        samfile = pysam.AlignmentFile(sf_bam, "rb", reference_filename=self.sf_reference)
        m_clip_pos = {}
        m_chrm_id_name = self._get_chrm_id_name(samfile)
        for algnmt in samfile.fetch(chrm):  ##fetch reads mapped to "chrm"
            ##here need to skip the secondary and supplementary alignments?
            # if algnmt.is_secondary or algnmt.is_supplementary:
            #     continue
            if algnmt.is_duplicate == True:  ##duplicate
                continue
            b_first = True
            if algnmt.is_read2 == True:
                b_first = False
            if algnmt.is_unmapped == True:  # unmapped
                continue
            l_cigar = algnmt.cigar
            if len(l_cigar) < 1:  # wrong alignment
                continue
            if len(l_cigar) == 1 and l_cigar[0][0] == 0:  ##fully mapped
                continue
            if algnmt.mapping_quality < global_values.MINIMUM_CLIP_MAPQ:#by default this is set to 12
                continue
####
            if algnmt.next_reference_id not in m_chrm_id_name:
                continue
            
            if len(l_cigar) > 2: # YW 2022/11/13: filter out multiple soft/hard clipped alignments
                num_clipped = 0
                for entry in l_cigar:
                    if entry[0] == 4 or entry[0] == 5:
                        num_clipped += 1
                if num_clipped > 1:
                    continue
                
            map_pos = algnmt.reference_start
            mate_chrm = '*'
            mate_pos = 0

            if algnmt.mate_is_unmapped == False and algnmt.next_reference_id>=0:
                mate_chrm = algnmt.next_reference_name
                mate_pos = algnmt.next_reference_start
            # print mate_chrm, mate_pos ################################################################################
            b_mate_in_rep_Alu, b_mate_in_rep_L1, b_mate_in_rep_SVA = False, False, False
            rep_start_pos_Alu, rep_start_pos_L1, rep_start_pos_SVA = 0, 0, 0
            #b_mate_in_rep, rep_start_pos = xannotation.is_within_repeat_region(mate_chrm, mate_pos)
            if b_se == False:
                b_mate_in_rep_Alu, rep_start_pos_Alu = xannotation_Alu.is_within_repeat_region_interval_tree(mate_chrm, mate_pos)
                b_mate_in_rep_L1, rep_start_pos_L1 = xannotation_L1.is_within_repeat_region_interval_tree(mate_chrm, mate_pos)
                b_mate_in_rep_SVA, rep_start_pos_SVA = xannotation_SVA.is_within_repeat_region_interval_tree(mate_chrm, mate_pos)
            
            # print b_mate_in_rep, rep_start_pos #######################################################################
            if l_cigar[0][0] == 4 or l_cigar[0][0] == 5:  # left clipped
                if map_pos not in m_clip_pos:
                    m_clip_pos[map_pos] = [1, 0, 0, 0, 0]
                else:
                    m_clip_pos[map_pos][0] += 1
                if b_mate_in_rep_Alu:
                    m_clip_pos[map_pos][2] += 1
                if b_mate_in_rep_L1:
                    m_clip_pos[map_pos][3] += 1
                if b_mate_in_rep_SVA:
                    m_clip_pos[map_pos][4] += 1

                if algnmt.is_supplementary or algnmt.is_secondary:  ###secondary and supplementary are not considered
                    continue

            if l_cigar[-1][0] == 4 or l_cigar[-1][0] == 5:  # right clipped
                ##calculate the exact clip position
                for (type, lenth) in l_cigar[:-1]:
                    if type == 4 or type == 5 or type == 1:  # (1 for insertion)
                        continue
                    else:
                        map_pos += lenth

                if map_pos not in m_clip_pos:
                    m_clip_pos[map_pos] = [0, 1, 0, 0, 0]
                else:
                    m_clip_pos[map_pos][1] += 1
                if b_mate_in_rep_Alu:
                    m_clip_pos[map_pos][2] += 1
                if b_mate_in_rep_L1:
                    m_clip_pos[map_pos][3] += 1
                if b_mate_in_rep_SVA:
                    m_clip_pos[map_pos][4] += 1

        sf_clip_pos = working_folder + chrm + global_values.CLIP_POS_SUFFIX
        with open(sf_clip_pos, "w") as fout_clip_pos:
            for pos in m_clip_pos:
                i_all_clip = m_clip_pos[pos][0] + m_clip_pos[pos][1]
                if i_all_clip >= i_clip_cutoff:  ####set a cutoff for the total number of clipped ones
                    fout_clip_pos.write("\t".join([str(pos), str(m_clip_pos[pos][0]), str(m_clip_pos[pos][1]),
                                                   str(m_clip_pos[pos][2]), str(m_clip_pos[pos][3]), str(m_clip_pos[pos][4])+"\n"]))
        samfile.close()
    
    
    # # YW 2021/10/07 to modularize collect_clip_info_and_parts_by_chrm, so we can fix errors with reads with left and right clipped sequences and left clip seq not passing the filter
    # def _filter_clipped_parts(self, len_clip_seq, clipped_seq, l_quality_score):
    #     # left clipped
    #     xpolyA=PolyA()
    #     if len_clip_seq>=global_values.MINIMUM_POLYA_CLIP and len_clip_seq<=global_values.BWA_REALIGN_CUTOFF:
    #         #check whether this the polyA side
    #         #Here have more strict requirement --> YW 2020/06/30: change this to be less stringent?
    #         # YW 2020/08/09 changing contain_enough_A_T to allow 25% non A-T
    #         min_AT_ratio=global_values.MIN_AT_RATIO
    #         if xpolyA.contain_enough_A_T(clipped_seq, len_clip_seq, min_AT_ratio)==False:
    #             return False
    #         # if xpolyA.contain_enough_A_T_old(clipped_seq, len_clip_seq)==False:
    #         #     continue
    #     elif len_clip_seq < global_values.BWA_REALIGN_CUTOFF:
    #         return False
    
    #     #require the medium quality score be higher than threshold
    #     if self._is_qualified_clip(l_quality_score)==False:
    #         return False
    #     return True
    
    def _write_clipped_parts(self, L_R, b_first, l_query_quality, query_name, mate_chrm, mate_pos, chrm, map_pos, f_clip_fq, clipped_seq):
        # left clipped
        clipped_qulity = self._cvt_to_Ascii_quality(l_query_quality)
        clipped_rname = global_values.SEPARATOR.join([query_name, mate_chrm, str(mate_pos), ""])
        if L_R == "L":
            s_tmp = "{0}{1}{2}{3}{4}{5}2".format(chrm, global_values.SEPARATOR, map_pos, global_values.SEPARATOR,
                                             global_values.FLAG_LEFT_CLIP, global_values.SEPARATOR)
            if b_first:
                s_tmp = "{0}{1}{2}{3}{4}{5}1".format(chrm, global_values.SEPARATOR, map_pos, global_values.SEPARATOR,
                                                 global_values.FLAG_LEFT_CLIP, global_values.SEPARATOR)
        else:
            s_tmp = "{0}{1}{2}{3}{4}{5}2".format(chrm, global_values.SEPARATOR, map_pos, global_values.SEPARATOR,
                                             global_values.FLAG_RIGHT_CLIP, global_values.SEPARATOR)
            if b_first:
                s_tmp = "{0}{1}{2}{3}{4}{5}1".format(chrm, global_values.SEPARATOR, map_pos, global_values.SEPARATOR,
                                                 global_values.FLAG_RIGHT_CLIP, global_values.SEPARATOR)
        clipped_rname += s_tmp

        f_clip_fq.write("".join(["@", clipped_rname, "\n"]))
        f_clip_fq.write("".join([clipped_seq, "\n+\n"]))
        f_clip_fq.write("".join([clipped_qulity, "\n"]))
        
        
    # YW 2021/09/29 added for predefined loci
    # YW 2021/10/07 fixed errors with reads with left and right clipped sequences (if left clipped seq shows up but with insufficient clip length, the right clipped seq would be skipped)
    def collect_clip_info_and_parts_by_chrm(self, record):
        chrm = record[0]
        pos_list = record[1]
        sf_bam = record[2]
        working_folder = record[3]
        sf_annotation_Alu = record[4]
        sf_annotation_L1 = record[5]
        sf_annotation_SVA = record[6]
        b_with_chr = record[7]
        i_clip_cutoff = int(record[8])
        b_se = record[9] #whether this is single end reads
        sf_clip_fq = working_folder + chrm + global_values.CLIP_FQ_SUFFIX  # this is to save the clipped part for re-alignment
        f_clip_fq = open(sf_clip_fq, "w")
####
        xannotation_Alu = XAnnotation(sf_annotation_Alu)
        xannotation_Alu.set_with_chr(b_with_chr)
        xannotation_L1 = XAnnotation(sf_annotation_L1)
        xannotation_L1.set_with_chr(b_with_chr)
        xannotation_SVA = XAnnotation(sf_annotation_SVA)
        xannotation_SVA.set_with_chr(b_with_chr)
        #xannotation.load_rmsk_annotation()
        i_min_copy_len=0
        i_boundary_extnd=0
        # if global_values.IS_CALL_SVA==True:
        #     i_boundary_extnd=global_values.SVA_ANNOTATION_EXTND
        xannotation_Alu.load_rmsk_annotation_with_extnd_with_lenth_cutoff(i_boundary_extnd, i_min_copy_len)
        xannotation_Alu.index_rmsk_annotation_interval_tree()
        xannotation_L1.load_rmsk_annotation_with_extnd_with_lenth_cutoff(i_boundary_extnd, i_min_copy_len)
        xannotation_L1.index_rmsk_annotation_interval_tree()
        xannotation_SVA.load_rmsk_annotation_with_extnd_with_lenth_cutoff(global_values.SVA_ANNOTATION_EXTND, i_min_copy_len)
        xannotation_SVA.index_rmsk_annotation_interval_tree()
        xpolyA=PolyA()
        
        
        samfile = pysam.AlignmentFile(sf_bam, "rb", reference_filename=self.sf_reference)
        m_clip_pos = {}
        m_chrm_id_name = self._get_chrm_id_name(samfile)
        for insertion_pos in pos_list:
            try: # 2021/11/29 added this try except to deal with input with both start and end positions
                int(insertion_pos)
                start_pos = insertion_pos - global_values.NEARBY_REGION
                end_pos = insertion_pos + global_values.NEARBY_REGION + 1
            except TypeError:
                start_pos, end_pos = insertion_pos
            if start_pos <= 0:
                start_pos = 1
            
            for algnmt in samfile.fetch(chrm, start_pos, end_pos):  ##fetch reads mapped to "chrm"
                ##here need to skip the secondary and supplementary alignments?
                # YW 2021/10/02 uncomment the first if statement below
                if algnmt.is_secondary or algnmt.is_supplementary:
                    continue
                if algnmt.is_duplicate == True:  ##duplicate
                    continue
                b_first = True
                if algnmt.is_read2 == True:
                    b_first = False
                if algnmt.is_unmapped == True:  # unmapped
                    continue
                l_cigar = algnmt.cigar
                if len(l_cigar) < 1:  # wrong alignment
                    continue
                if len(l_cigar) == 1 and l_cigar[0][0] == 0:  ##fully mapped
                    continue
                if algnmt.mapping_quality < global_values.MINIMUM_CLIP_MAPQ:#by default this is set to 12
                    continue
    ####
                if algnmt.next_reference_id not in m_chrm_id_name:
                    continue
                
                map_pos = algnmt.reference_start
                mate_chrm = '*'
                mate_pos = 0
    
                if algnmt.mate_is_unmapped == False and algnmt.next_reference_id>=0:
                    mate_chrm = algnmt.next_reference_name
                    mate_pos = algnmt.next_reference_start
                
                # YW 2021/10/02 added from collect_clipped_parts_by_chrm
                query_name = algnmt.query_name
                query_seq = algnmt.query_sequence
                query_quality = algnmt.query_qualities
                # print mate_chrm, mate_pos ################################################################################
                b_mate_in_rep_Alu, b_mate_in_rep_L1, b_mate_in_rep_SVA = False, False, False
                rep_start_pos_Alu, rep_start_pos_L1, rep_start_pos_SVA = 0, 0, 0
                #b_mate_in_rep, rep_start_pos = xannotation.is_within_repeat_region(mate_chrm, mate_pos)
                if b_se == False:
                    b_mate_in_rep_Alu, rep_start_pos_Alu = xannotation_Alu.is_within_repeat_region_interval_tree(mate_chrm, mate_pos)
                    b_mate_in_rep_L1, rep_start_pos_L1 = xannotation_L1.is_within_repeat_region_interval_tree(mate_chrm, mate_pos)
                    b_mate_in_rep_SVA, rep_start_pos_SVA = xannotation_SVA.is_within_repeat_region_interval_tree(mate_chrm, mate_pos)
    
                # print b_mate_in_rep, rep_start_pos #######################################################################
                if l_cigar[0][0] == 4 or l_cigar[0][0] == 5:  # left clipped
                    if map_pos not in m_clip_pos:
                        m_clip_pos[map_pos] = [1, 0, 0, 0, 0]
                    else:
                        m_clip_pos[map_pos][0] += 1
                    if b_mate_in_rep_Alu:
                        m_clip_pos[map_pos][2] += 1
                    if b_mate_in_rep_L1:
                        m_clip_pos[map_pos][3] += 1
                    if b_mate_in_rep_SVA:
                        m_clip_pos[map_pos][4] += 1
                    
                    # YW 2021/10/02 commented out below
                    # if algnmt.is_supplementary or algnmt.is_secondary:  ###secondary and supplementary are not considered
                    #     continue
                    
                    #######################################################################
                    # YW 2021/10/02 copied below from collect_clipped_parts_by_chrm
                    clipped_seq = query_seq[:l_cigar[0][1]]
                    len_clip_seq=len(clipped_seq)
                    for ch in query_seq:
                        if ch=="N":
                            len_clip_seq -= 1
                    l_query_quality = query_quality[:l_cigar[0][1]]
                    pass_filter = self._filter_clipped_parts(len_clip_seq, clipped_seq, l_query_quality)
                    if pass_filter:
                        self._write_clipped_parts("L", b_first, l_query_quality, query_name, mate_chrm, mate_pos, chrm, map_pos, f_clip_fq, clipped_seq)
                    #######################################################################
    
                if l_cigar[-1][0] == 4 or l_cigar[-1][0] == 5:  # right clipped
                    ##calculate the exact clip position
                    for (type, lenth) in l_cigar[:-1]:
                        if type == 4 or type == 5 or type == 1:  # (1 for insertion)
                            continue
                        else:
                            map_pos += lenth
    
                    if map_pos not in m_clip_pos:
                        m_clip_pos[map_pos] = [0, 1, 0, 0, 0]
                    else:
                        m_clip_pos[map_pos][1] += 1
                    if b_mate_in_rep_Alu:
                        m_clip_pos[map_pos][2] += 1
                    if b_mate_in_rep_L1:
                        m_clip_pos[map_pos][3] += 1
                    if b_mate_in_rep_SVA:
                        m_clip_pos[map_pos][4] += 1
                    
                    #######################################################################
                    # YW 2021/10/02 copied below from collect_clipped_parts_by_chrm
                    start_pos = -1 * l_cigar[-1][1]
                    clipped_seq = query_seq[start_pos:]
                    len_clip_seq = len(clipped_seq)
                    for ch in query_seq:
                        if ch=="N":
                            len_clip_seq -= 1
                    
                    l_query_quality = query_quality[start_pos:]
                    pass_filter = self._filter_clipped_parts(len_clip_seq, clipped_seq, query_quality[start_pos:])
                    if pass_filter:
                        self._write_clipped_parts("R", b_first, l_query_quality, query_name, mate_chrm, mate_pos, chrm, map_pos, f_clip_fq, clipped_seq)
                    ##################################################################

        sf_clip_pos = working_folder + chrm + global_values.CLIP_POS_SUFFIX
        with open(sf_clip_pos, "w") as fout_clip_pos:
            for pos in m_clip_pos:
                i_all_clip = m_clip_pos[pos][0] + m_clip_pos[pos][1]
                if i_all_clip >= i_clip_cutoff:  ####set a cutoff for the total number of clipped ones
                    fout_clip_pos.write("\t".join([str(pos), str(m_clip_pos[pos][0]), str(m_clip_pos[pos][1]),
                                                   str(m_clip_pos[pos][2]), str(m_clip_pos[pos][3]), str(m_clip_pos[pos][4])+"\n"]))
        samfile.close()
        # YW 2021/10/02 copied below from collect_clipped_parts_by_chrm
        f_clip_fq.close()
####

    ####This function return:
    ####1. dictionary of clip position, ##in format {chrm: {map_pos: (left_cnt, right_cnt)}}
    def collect_clip_positions(self, sf_annotation_Alu, sf_annotation_L1, sf_annotation_SVA, i_clip_cutoff, b_se, sf_pub_folder):
        bam_info = BamInfo(self.sf_bam, self.sf_reference)
        b_with_chr = bam_info.is_chrm_contain_chr()
        samfile = pysam.AlignmentFile(self.sf_bam, "rb", reference_filename=self.sf_reference)
        references = samfile.references
        xchrom=XChromosome()
        l_chrm_records = []
        for chrm in references:
            if xchrom.is_decoy_contig_chrms(chrm) == True:  ###decoy sequnces and contigs are not considered
                continue
            l_chrm_records.append((chrm, self.sf_bam, self.working_folder,
                                   sf_annotation_Alu, sf_annotation_L1, sf_annotation_SVA,
                                   b_with_chr, i_clip_cutoff, b_se)) # YW 2021/03/18 update sf_annotation
        samfile.close()

        pool = Pool(self.n_jobs)
        pool.map(unwrap_self_cnt_clip_pos, list(zip([self] * len(l_chrm_records), l_chrm_records)), 1)
        pool.close()
        pool.join()

        #soft_link clip pos
        for rcd in l_chrm_records:
            sf_clip_pos = self.working_folder + rcd[0] + global_values.CLIP_POS_SUFFIX
            sf_pub_pos=sf_pub_folder + rcd[0] + global_values.CLIP_POS_SUFFIX
            if os.path.islink(sf_pub_pos)==True or os.path.isfile(sf_pub_pos)==True:
                os.remove(sf_pub_pos)
            cmd="ln -s {0} {1}".format(sf_clip_pos, sf_pub_folder)
            Popen(cmd, shell=True, stdout=PIPE).communicate()
    
    
    # YW 2021/09/29 new function, only collect info and clipped parts for predefined loci
    def collect_clip_info_and_parts(self, locus_dict, sf_annotation_Alu, sf_annotation_L1, sf_annotation_SVA, i_clip_cutoff, b_se, sf_pub_folder, sf_all_clip_fq_ori):
        bam_info = BamInfo(self.sf_bam, self.sf_reference)
        b_with_chr = bam_info.is_chrm_contain_chr()
        samfile = pysam.AlignmentFile(self.sf_bam, "rb", reference_filename=self.sf_reference)
        references = samfile.references
        xchrom=XChromosome()
        l_chrm_records = []
        for chrm in locus_dict:
            if xchrom.is_decoy_contig_chrms(chrm) == True:  ###decoy sequnces and contigs are not considered
                continue
            l_chrm_records.append((chrm, locus_dict[chrm], self.sf_bam, self.working_folder,
                                   sf_annotation_Alu, sf_annotation_L1, sf_annotation_SVA,
                                   b_with_chr, i_clip_cutoff, b_se)) # YW 2021/03/18 update sf_annotation
        samfile.close()

        pool = Pool(self.n_jobs)
        pool.map(unwrap_self_collect_clip_info_and_parts, list(zip([self] * len(l_chrm_records), l_chrm_records)), 1)
        pool.close()
        pool.join()

        #soft_link clip pos
        for rcd in l_chrm_records:
            sf_clip_pos = self.working_folder + rcd[0] + global_values.CLIP_POS_SUFFIX
            sf_pub_pos=sf_pub_folder + rcd[0] + global_values.CLIP_POS_SUFFIX
            if os.path.islink(sf_pub_pos)==True or os.path.isfile(sf_pub_pos)==True:
                os.remove(sf_pub_pos)
            cmd="ln -s {0} {1}".format(sf_clip_pos, sf_pub_folder)
            Popen(cmd, shell=True, stdout=PIPE).communicate()
            
        with open(sf_all_clip_fq_ori, "w") as fout_all:
            for chrm in references:
                sf_clip_fq = self.working_folder + chrm + global_values.CLIP_FQ_SUFFIX
                if os.path.isfile(sf_clip_fq) == False:
                    continue
                with open(sf_clip_fq) as fin_clip:
                    for line in fin_clip:
                        fout_all.write(line)
                os.remove(sf_clip_fq)#clean the temporary file

    ####given specific chrm, get all the related clipped reads
    #Here besides the long clipped parts, we keep the very short clipped parts that have dominant A or T
    def collect_clipped_parts_by_chrm(self, record):
        chrm = record[0]
        sf_bam = record[1]
        working_folder = record[2]

        # first load in the positions by chromosome
        sf_clip_pos = working_folder + chrm + global_values.CLIP_POS_SUFFIX
        m_pos = {}
        if os.path.exists(sf_clip_pos) == False:
            return
        with open(sf_clip_pos) as fin_clip_pos:
            for line in fin_clip_pos:
                fields = line.split()
                pos = int(fields[0])
                m_pos[pos] = 1
        # second, load the reads, and write the related clipped part into file
        sf_clip_fq = working_folder + chrm + global_values.CLIP_FQ_SUFFIX  # this is to save the clipped part for re-alignment
        f_clip_fq = open(sf_clip_fq, "w")
        xpolyA=PolyA()
        samfile = pysam.AlignmentFile(sf_bam, "rb", reference_filename=self.sf_reference)
        m_chrm_id=self._get_chrm_id_name(samfile)
        for algnmt in samfile.fetch(chrm):  ##fetch reads mapped to "chrm"
            if algnmt.is_secondary or algnmt.is_supplementary:
                continue
            if algnmt.is_duplicate == True:  ##duplicate
                continue
            b_first = True
            if algnmt.is_read2 == True:
                b_first = False
            if algnmt.is_unmapped == True:  # unmapped
                continue
            l_cigar = algnmt.cigar
            if len(l_cigar) < 1:  # wrong alignment
                continue
            if len(l_cigar) == 1 and l_cigar[0][0] == 0:  ##fully mapped
                continue

            query_name = algnmt.query_name
            query_seq = algnmt.query_sequence
            query_quality = algnmt.query_qualities  ##this is different from the one saved in the fastq/sam, no offset 33 to subtract
            map_pos = algnmt.reference_start
            mate_chrm = '*'
            mate_pos = 0
            if algnmt.next_reference_id not in m_chrm_id:
                continue
            if algnmt.mate_is_unmapped == False and algnmt.next_reference_id>=0:
                mate_chrm = algnmt.next_reference_name
                mate_pos = algnmt.next_reference_start
            
            if len(l_cigar) > 2: # 2022/11/13 YW: filter out multiple soft clipped alignments
                num_clipped = 0
                for entry in l_cigar:
                    if entry[0] == 4 or entry[0] == 5:
                        num_clipped += 1
                if num_clipped > 1:
                    continue
                
            if l_cigar[0][0] == 4:  # left soft clipped
                if map_pos in m_pos:
                    clipped_seq = query_seq[:l_cigar[0][1]]
                    # YW 2020/08/09 added this to account for mask 2
                    len_clip_seq=len(clipped_seq)
                    for ch in algnmt.query_sequence:
                        if ch=="N":
                            len_clip_seq -= 1
                    if len_clip_seq>=global_values.MINIMUM_POLYA_CLIP and len_clip_seq<=global_values.BWA_REALIGN_CUTOFF:
                        #check whether this the polyA side
                        #Here have more strict requirement --> YW 2020/06/30: change this to be less stringent?
                        # YW 2020/08/09 changing contain_enough_A_T to allow 25% non A-T
                        min_AT_ratio=global_values.MIN_AT_RATIO
                        if xpolyA.contain_enough_A_T(clipped_seq, len_clip_seq, min_AT_ratio)==False:
                            continue
                        # if xpolyA.contain_enough_A_T_old(clipped_seq, len_clip_seq)==False:
                        #     continue
                    elif len_clip_seq < global_values.BWA_REALIGN_CUTOFF:
                        continue

                    #require the medium quality score should higher than threshold
                    l_quality_score=query_quality[:l_cigar[0][1]]
                    if self._is_qualified_clip(l_quality_score)==False:
                        continue
                    clipped_qulity = self._cvt_to_Ascii_quality(query_quality[:l_cigar[0][1]])
                    clipped_rname = global_values.SEPARATOR.join([query_name, mate_chrm, str(mate_pos), ""])
                    s_tmp = "{0}{1}{2}{3}{4}{5}2".format(chrm, global_values.SEPARATOR, map_pos, global_values.SEPARATOR,
                                                         global_values.FLAG_LEFT_CLIP, global_values.SEPARATOR)
                    clipped_rname += s_tmp
                    if b_first:
                        clipped_rname = global_values.SEPARATOR.join([query_name, mate_chrm, str(mate_pos), ""])
                        s_tmp = "{0}{1}{2}{3}{4}{5}1".format(chrm, global_values.SEPARATOR, map_pos, global_values.SEPARATOR,
                                                             global_values.FLAG_LEFT_CLIP, global_values.SEPARATOR)
                        clipped_rname += s_tmp

                    f_clip_fq.write("".join(["@", clipped_rname, "\n"]))
                    f_clip_fq.write("".join([clipped_seq, "\n+\n"]))
                    f_clip_fq.write("".join([clipped_qulity, "\n"]))

            if l_cigar[-1][0] == 4:  #right clipped
                ##calculate the exact clip position
                for (type, lenth) in l_cigar[:-1]:
                    if type == 4 or type == 5 or type == 1:  # (1 for insertion)
                        continue
                    else:
                        map_pos += lenth

                if map_pos in m_pos:  # soft-clip
                    clipped_rname = global_values.SEPARATOR.join([query_name, mate_chrm, str(mate_pos), ""])
                    s_tmp = "{0}{1}{2}{3}{4}{5}2".format(chrm, global_values.SEPARATOR, map_pos, global_values.SEPARATOR,
                                                         global_values.FLAG_RIGHT_CLIP, global_values.SEPARATOR)
                    clipped_rname += s_tmp
                    start_pos = -1 * l_cigar[-1][1]
                    clipped_seq = query_seq[start_pos:]
                    len_clip_seq = len(clipped_seq)
                    for ch in algnmt.query_sequence:
                        if ch=="N":
                            len_clip_seq -= 1
                    if len_clip_seq >= global_values.MINIMUM_POLYA_CLIP and len_clip_seq <= global_values.BWA_REALIGN_CUTOFF:
                        # YW 2020/08/09 changing contain_enough_A_T to allow 25% nonA-T
                        min_AT_ratio=global_values.MIN_AT_RATIO
                        if xpolyA.contain_enough_A_T(clipped_seq, len_clip_seq, min_AT_ratio)==False:
                            continue
                        # if xpolyA.contain_enough_A_T_old(clipped_seq, len_clip_seq)==False:
                        #     continue
                    elif len_clip_seq < global_values.BWA_REALIGN_CUTOFF:
                        continue

                    l_quality_score = query_quality[start_pos:]
                    if self._is_qualified_clip(l_quality_score) == False:
                        continue

                    clipped_qulity = self._cvt_to_Ascii_quality(query_quality[start_pos:])
                    if b_first:
                        clipped_rname = global_values.SEPARATOR.join([query_name, mate_chrm, str(mate_pos), ""])
                        s_tmp = "{0}{1}{2}{3}{4}{5}1".format(chrm, global_values.SEPARATOR, map_pos, global_values.SEPARATOR,
                                                             global_values.FLAG_RIGHT_CLIP, global_values.SEPARATOR)
                        clipped_rname += s_tmp

                    f_clip_fq.write("".join(["@", clipped_rname, "\n"]))
                    f_clip_fq.write("".join([clipped_seq, "\n+\n"]))
                    f_clip_fq.write("".join([clipped_qulity, "\n"]))
        samfile.close()
        f_clip_fq.close()

    ####This function: YW 2020/07/19 updated the comments
    ####1. given TE/tmp/clip/chr.clip_pos, format: pos, left_cnt, right_cnt, num_mate_in_rep (now add Alu, L1, SVA separately)
    ####2. do some filtering
    ####3. write the clipped part into fastq files
    def collect_clipped_parts(self, sf_all_clip_fq):
        samfile = pysam.AlignmentFile(self.sf_bam, "rb", reference_filename=self.sf_reference)
        references = samfile.references
        l_chrm_records = []
        xchrom = XChromosome()
        for chrm in references:
            if xchrom.is_decoy_contig_chrms(chrm) == True:  ###Here filter out those aligned to decoy sequences
                #print(("Skip chromosome {0}".format(chrm)))  ##YW 2020/07/21 disable printing
                continue
            l_chrm_records.append((chrm, self.sf_bam, self.working_folder))
        samfile.close()

        pool = Pool(self.n_jobs)
        pool.map(unwrap_self_collect_clip_parts, list(zip([self] * len(l_chrm_records), l_chrm_records)), 1)
        pool.close()
        pool.join()

        ##merge the clipped reads
        with open(sf_all_clip_fq, "w") as fout_all:
            for chrm in references:
                sf_clip_fq = self.working_folder + chrm + global_values.CLIP_FQ_SUFFIX
                if os.path.isfile(sf_clip_fq) == False:
                    continue
                with open(sf_clip_fq) as fin_clip:
                    for line in fin_clip:
                        fout_all.write(line)
                os.remove(sf_clip_fq)#clean the temporary file


    ####
    def run_cnt_clip_part_aligned_to_rep_by_chrm(self, record):
        ref_chrm = record[0]
        sf_sam_Alu = record[1]
        sf_sam_L1 = record[2]
        sf_sam_SVA = record[3]
        working_folder = record[4]

        ####clip positions for specific chrm
        # sf_clip_pos = working_folder + ref_chrm + global_values.CLIP_POS_SUFFIX
        # if os.path.isfile(sf_clip_pos) == False:
        #     print "Error: Position file for chrom {0} doesn't exist!!!!".format(ref_chrm)
        #     return
        ###clip positions for specific chrm with extra #left_clip #right_clip (mapped parts)
        sf_out_clip_pos = working_folder + ref_chrm + global_values.CLIP_RE_ALIGN_POS_SUFFIX

        m_sites_chrm = {}
        def run_cnt_clip_part_aligned_to_rep_by_chrm_helper(sf_sam, rep_type): # YW 2021/03/18 added to avoid duplicating code
            nonlocal m_sites_chrm # so we can change m_sites_chrm
            samfile = pysam.AlignmentFile(sf_sam, "r")
            for algnmt in samfile.fetch():  ##fetch reads mapped to "chrm"
                ##here need to skip the secondary and supplementary alignments?
                if algnmt.is_secondary or algnmt.is_supplementary:
                    continue
                # if algnmt.is_duplicate == True:  ##duplicate
                #     continue
                if algnmt.is_unmapped == True:
                    continue
    
                qname = algnmt.query_name
                qname_fields = qname.split(global_values.SEPARATOR)
                # chrm, map_pos, global_values.FLAG_LEFT_CLIP, first_read
                ori_chrm = qname_fields[-4]  #############check reverse-complementary consistent here ??????????????
                ori_mpos = int(qname_fields[-3])
    
                if ori_chrm != ref_chrm:  ###not the interesting chromosome
                    continue
    
                l_cigar = algnmt.cigar
                if len(l_cigar) < 1:  # wrong alignment
                    continue
                if len(l_cigar) > 2:
                    ####check the cigar
                    ###if both clipped, and the clipped part is large, then skip
                    b_left_clip = False
                    i_left_clip_len = 0
                    if l_cigar[0][0] == 4 or l_cigar[0][0] == 5:  # left clipped
                        b_left_clip = True
                        i_left_clip_len = l_cigar[0][1]
                    b_right_clip = False
                    i_right_clip_len = 0
                    if l_cigar[-1][0] == 4 or l_cigar[-1][0] == 5:  # right clipped
                        b_right_clip = True
                        i_right_clip_len = l_cigar[-1][1]
                    
                    if b_left_clip == True and b_right_clip == True:
                        if (i_left_clip_len > global_values.MAX_CLIP_CLIP_LEN) and (i_right_clip_len > global_values.MAX_CLIP_CLIP_LEN):
                            continue
    
                ####for the alignment (of the clipped read), if the mapped part is smaller than the clipped part,
                ####then skip
                n_total = 0
                n_map = 0
                for (type, lenth) in l_cigar:
                    if type == 0:
                        n_map += lenth
                    if type != 2:  # deletion is not added to the total length data
                        n_total += lenth
                
                # YW 2020/08/09 added this to account for mask 2
                for ch in algnmt.query_sequence:
                    if ch=="N":
                        n_total -= 1
        
                if n_map < (n_total * 3 / 4):  ########################require at least 3/4 of the seq is mapped !!!!!!!!!!!
                    continue
    
                b_left = True
                if qname_fields[-2] == global_values.FLAG_RIGHT_CLIP:
                    b_left = False
                
                # merge counts from Alu, L1, SVA sam files
                if (ori_mpos in m_sites_chrm) == False:
                    m_sites_chrm[ori_mpos] = []
                    if rep_type == "Alu":
                        if b_left == True:
                            m_sites_chrm[ori_mpos].extend([1, 0, 0, 0, 0, 0])
                        else:
                            m_sites_chrm[ori_mpos].extend([0, 1, 0, 0, 0, 0])
                    elif rep_type == "L1":
                        if b_left == True:
                            m_sites_chrm[ori_mpos].extend([0, 0, 1, 0, 0, 0])
                        else:
                            m_sites_chrm[ori_mpos].extend([0, 0, 0, 1, 0, 0])
                    elif rep_type == "SVA":
                        if b_left == True:
                            m_sites_chrm[ori_mpos].extend([0, 0, 0, 0, 1, 0])
                        else:
                            m_sites_chrm[ori_mpos].extend([0, 0, 0, 0, 0, 1])
                else:
                    if rep_type == "Alu":
                        if b_left == True:
                            m_sites_chrm[ori_mpos][0] += 1
                        else:
                            m_sites_chrm[ori_mpos][1] += 1
                    elif rep_type == "L1":
                        if b_left == True:
                            m_sites_chrm[ori_mpos][2] += 1
                        else:
                            m_sites_chrm[ori_mpos][3] += 1
                    elif rep_type == "SVA":
                        if b_left == True:
                            m_sites_chrm[ori_mpos][4] += 1
                        else:
                            m_sites_chrm[ori_mpos][5] += 1      
            samfile.close()
        run_cnt_clip_part_aligned_to_rep_by_chrm_helper(sf_sam_Alu, "Alu")
        run_cnt_clip_part_aligned_to_rep_by_chrm_helper(sf_sam_L1, "L1")
        run_cnt_clip_part_aligned_to_rep_by_chrm_helper(sf_sam_SVA, "SVA")
        
        
        with open(sf_out_clip_pos, "w") as fout_clip_pos:
            for pos in m_sites_chrm:
                # lth = len(m_sites_chrm[pos])
                fout_clip_pos.write(str(pos) + "\t")
                fout_clip_pos.write("\t".join([str(i) for i in m_sites_chrm[pos]])+"\n")
                # for i in range(lth):
                #     fout_clip_pos.write(str(m_sites_chrm[pos][i]) + "\t")
                # fout_clip_pos.write("\n")

    ####Input:1. the re-aligned clipped parts
    def cnt_clip_part_aligned_to_rep(self, sf_clip_sam_Alu, sf_clip_sam_L1, sf_clip_sam_SVA):
        samfile = pysam.AlignmentFile(self.sf_bam, "rb", reference_filename=self.sf_reference)
        references = samfile.references
        l_chrm_records = []
        #here filter out those uninterested contigs
        xchrom=XChromosome()
        for chrm in references:
            if xchrom.is_decoy_contig_chrms(chrm) == True:  ###filter out decoy and other contigs
                continue
            l_chrm_records.append((chrm, sf_clip_sam_Alu, sf_clip_sam_L1, sf_clip_sam_SVA, self.working_folder))
        samfile.close()

        pool = Pool(self.n_jobs)
        pool.map(unwrap_self_cnt_clip_parts, list(zip([self] * len(l_chrm_records), l_chrm_records)), 1)
        pool.close()
        pool.join()


    def merge_clip_positions(self, sf_pclip, sf_out):
        samfile = pysam.AlignmentFile(self.sf_bam, "rb", reference_filename=self.sf_reference)
        references = samfile.references

        with open(sf_out, "w") as fout_clip_pos:
            for chrm in references:
                ##first, load in the whole file into dict
                m_realign_pos = {}
                sf_chrm_re_align_clip_pos = self.working_folder + chrm + global_values.CLIP_RE_ALIGN_POS_SUFFIX
                if os.path.isfile(sf_chrm_re_align_clip_pos) == True:
                    with open(sf_chrm_re_align_clip_pos) as fin_realign_clip_chrm:
                        for line in fin_realign_clip_chrm:
                            fields = line.split()
                            m_realign_pos[int(fields[0])] = "\t".join(fields[1:])

                sf_clip_pos = sf_pclip + chrm + global_values.CLIP_POS_SUFFIX
                if os.path.isfile(sf_clip_pos) == False:
                    continue
                with open(sf_clip_pos) as fin_clip_pos:
                    for line in fin_clip_pos:
                        fields = line.split()
                        pos = int(fields[0])
                        if pos in m_realign_pos:
                            fout_clip_pos.write(chrm + "\t")
                            fout_clip_pos.write(line.rstrip() + "\t")
                            fout_clip_pos.write(m_realign_pos[pos] + "\n")
                        else:
                            fout_clip_pos.write(chrm + "\t")
                            fout_clip_pos.write(line.rstrip() + "\t")
                            fout_clip_pos.write("0\t0\n")
                #os.remove(sf_clip_pos)
        samfile.close()
####
    ####
    def merge_clip_positions_with_cutoff(self, cutoff_left_clip, cutoff_right_clip, max_cov_cutoff, sf_pclip, sf_out):
        samfile = pysam.AlignmentFile(self.sf_bam, "rb", reference_filename=self.sf_reference)
        references = samfile.references

        with open(sf_out, "w") as fout_clip_pos:
            for chrm in references:
                ##first, load in the whole file into dict
                m_realign_pos = {}
                sf_chrm_re_align_clip_pos = self.working_folder + chrm + global_values.CLIP_RE_ALIGN_POS_SUFFIX
                if os.path.isfile(sf_chrm_re_align_clip_pos) == True:
                    with open(sf_chrm_re_align_clip_pos) as fin_realign_clip_chrm:
                        for line in fin_realign_clip_chrm:
                            fields = line.split()
                            m_realign_pos[int(fields[0])] = "\t".join(fields[1:])
                    os.remove(sf_chrm_re_align_clip_pos) # YW 2021/04/07 remove unused files

                sf_clip_pos = sf_pclip + chrm + global_values.CLIP_POS_SUFFIX
                if os.path.isfile(sf_clip_pos) == False:
                    continue
                with open(sf_clip_pos) as fin_clip_pos:
                    for line in fin_clip_pos:
                        fields = line.split()
                        pos = int(fields[0])
                        i_left_clip = int(fields[1])
                        i_right_clip = int(fields[2])
                        if i_left_clip < cutoff_left_clip and i_right_clip < cutoff_right_clip:
                            continue
                        if (i_left_clip+i_right_clip) > max_cov_cutoff:
                            continue
                        if pos in m_realign_pos:
                            fout_clip_pos.write(chrm + "\t")
                            fout_clip_pos.write(line.rstrip() + "\t")
                            fout_clip_pos.write(m_realign_pos[pos] + "\n")
                        else:
                            fout_clip_pos.write(chrm + "\t")
                            fout_clip_pos.write(line.rstrip() + "\t")
                            fout_clip_pos.write("\t".join([str(0)]*6) + "\n")
                os.remove(sf_clip_pos) # YW uncommented 2021/04/07
        samfile.close()

