##11/22/2017
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu
'''
# defined the following classes: TE_Multi_Locator, TELocator
# YW 2020/08/03 github update:
# filter_candidate_sites_by_discordant_pairs_multi_alignmts function, filter_candidate_sites_by_discordant_pairs_non_barcode function, run_filter_by_discordant_pair_by_chrom_non_barcode
'''
import os
import sys
import pysam
from subprocess import *
from multiprocessing import Pool
from time import sleep # YW 2021/05/26 added for parallel realignment of clip/disc reads (to wait for cns realignment to finish for L1 and SVA before proceeding)

from clip_read import ClipReadInfo
from x_annotation import *
from x_alignments import *
from x_intermediate_sites import *
from bwa_align import *
import global_values
from cmd_runner import *
from cns_parallel import *  # YW 2021/05/26 added for parallel realignment of clip/disc reads


# YW 2021/04/21 to write lists into final output
def conv_list_to_str(lst):
	if len(lst) == 0:
		string = "[]"
	string = "["+",".join(list(map(str, lst)))+"]"
	return string


class TE_Multi_Locator():
	def __init__(self, sf_list, s_working_folder, n_jobs, sf_ref):
		self.sf_list = sf_list
		self.working_folder = s_working_folder
		self.n_jobs = int(n_jobs)
		self.sf_ref=sf_ref ##reference genome
		self.cmd_runner = CMD_RUNNER()

	def get_clip_part_realignment_for_list(self, sf_candidate_sites):
		m_sites={}
		with open(sf_candidate_sites) as fin_sites:
			for line in fin_sites:
				fields=line.split()
				chrm=fields[0]
				pos=int(fields[1])
				if chrm not in m_sites:
					m_sites[chrm]={}
				m_sites[chrm][pos]=line ###here save the chrom and position, also the other informations
		
		#for each bam in the list
		with open(self.sf_list) as fin_bam_list:
			for line in fin_bam_list:  ###for each bam file
				fields=line.split()
				sf_ori_bam = fields[0]
				if len(sf_ori_bam) <= 1:
					continue

	def collect_all_clipped_from_multiple_alignmts(self, sf_annotation, b_se, s_clip_wfolder, wfolder_pub_clip):
		with open(self.sf_list) as fin_bam_list:
			for line in fin_bam_list:  ###for each bam file
				fields = line.split()
				sf_ori_bam = fields[0]
				if len(sf_ori_bam) <= 1:
					continue
				
				caller = TELocator(sf_ori_bam, sf_ori_bam, self.working_folder, self.n_jobs, self.sf_ref)
				caller.collect_all_clipped_reads_only(sf_annotation, b_se, s_clip_wfolder, wfolder_pub_clip)
####

	####
	# YW 2020/08/01 added b_mosaic in input (github update)
	# YW 2021/05/19 took out b_mosaic
	def call_TEI_candidate_sites_from_multiple_alignmts(self, sf_annotation_Alu, sf_annotation_L1, sf_annotation_SVA,
							    sf_rep_cns_Alu, sf_rep_cns_L1, sf_rep_cns_SVA,
							    sf_rep_Alu, sf_rep_L1, sf_rep_SVA,
							    b_se, cutoff_left_clip,
							    cutoff_right_clip, cutoff_clip_mate_in_rep, cutoff_clip_mate_in_cns,
							    sf_clip_folder, b_force, max_cov, sf_out):
		cnt = 0
		s_sample_bam = ""
		b_set = False
		with open(self.sf_list) as fin_bam_list:
			i_idx_bam=0#indicates which bam this is
			for line in fin_bam_list:  ###for each bam file
				fields=line.split()
				if len(fields)<2:
					print("Error input bam file {0}!!!".format(line.rstrip()))
					continue
				sf_ori_bam = fields[0]
				s_read_type=fields[1].rstrip()
				print("Input bam {0} is sequenced from {1} platform!".format(sf_ori_bam, s_read_type))
				if global_values.X10 == s_read_type:###for 10X, we use a larger cutoff, as more clipped reads
					print("10X bam! Set the initial cutoff as {0}".format(global_values.INITIAL_MIN_CLIP_CUTOFF_X10))
					global_values.set_initial_min_clip_cutoff(global_values.INITIAL_MIN_CLIP_CUTOFF_X10)
				else:
					# YW 2020/08/09 want adjustment to default INITIAL_MIN_CLIP_CUTOFF_ILLUMINA regardless of b_mosaic
					# if cutoff_left_clip<=2 and b_mosaic==True :#
					if cutoff_left_clip<=2:
						print("Clip cutoff is small (<=2) , we are using 1 for initial cutoff")
						global_values.set_initial_min_clip_cutoff(1)#for low coverage data, set this to 1
					else:
						global_values.set_initial_min_clip_cutoff(global_values.INITIAL_MIN_CLIP_CUTOFF_ILLUMINA)
#
				if len(sf_ori_bam) <= 1:
					continue
				if b_set == False:
					s_sample_bam = sf_ori_bam
					b_set = True
				
				b_cutoff = True
				cutoff_hit_rep_copy=global_values.INITIAL_MIN_CLIP_CUTOFF
				
				# view the barcode bam as normal illumina bam
				# for each alignment, has one output
				sf_out_tmp = self.working_folder + global_values.CLIP_TMP + '{0}'.format(cnt)
				cnt += 1
				
				caller = TELocator(sf_ori_bam, sf_ori_bam, self.working_folder, self.n_jobs, self.sf_ref)
				# s_working_folder + global_values.CLIP_FOLDER + "/"+sf_bam_name + global_values.CLIP_FQ_SUFFIX
				sf_new_pub=""
				if len(sf_clip_folder)==0 or sf_clip_folder==None:
					print("public folder is null!!!!")
					continue
				if sf_clip_folder[-1]=="/":
					sf_new_pub=sf_clip_folder+"{0}/".format(i_idx_bam)
				else:
					sf_new_pub = sf_clip_folder + "/{0}/".format(i_idx_bam)
				caller.call_TEI_candidate_sites_from_clip_reads_v2(sf_annotation_Alu, sf_annotation_L1, sf_annotation_SVA,
										   sf_rep_cns_Alu, sf_rep_cns_L1, sf_rep_cns_SVA,
										   sf_rep_Alu, sf_rep_L1, sf_rep_SVA,
										   b_se, cutoff_hit_rep_copy, cutoff_hit_rep_copy, b_cutoff,
										   sf_new_pub, i_idx_bam, b_force, max_cov, sf_out_tmp)
				i_idx_bam+=1
####
		# get all the chromsomes names
		bam_info = BamInfo(s_sample_bam, self.sf_ref)
		b_with_chr = bam_info.is_chrm_contain_chr()
		m_chrms = bam_info.get_all_reference_names()
		
		xfilter = XIntermediateSites()
		xchrom=XChromosome()
		sf_out_merged = sf_out + "_tmp"
		with open(sf_out_merged, "w") as fout_sites_merged, open(sf_out, "w") as fout_sites:
			for chrm in m_chrms:  # write out chrm by chrm to save memory
				if xchrom.is_decoy_contig_chrms(chrm) == True:  ###filter out decoy and other contigs
					continue
				m_sites_chrm = {}
				for i in range(cnt):
					sf_tmp = self.working_folder + global_values.CLIP_TMP + "{0}".format(i)
					if os.path.isfile(sf_tmp) == False:
						print("Errors happen, file {0} doens't exist!".format(sf_tmp))
						continue
					with open(sf_tmp) as fin_tmp:
						for line in fin_tmp:
							fields = line.split()
							tmp_chrm = bam_info.process_chrm_name(fields[0], b_with_chr)
							if tmp_chrm != chrm:
								continue
							pos = int(fields[1])
							
							if pos not in m_sites_chrm:
								m_sites_chrm[pos] = []
								for value in fields[2:]:
									m_sites_chrm[pos].append(int(value))
							else:
								i_value = 0
								for value in fields[2:]:
									###sum (left-realign, right-realign, mate_in_rep)
									m_sites_chrm[pos][i_value] += int(value)
									i_value += 1

				for pos in m_sites_chrm:
					fout_sites_merged.write(chrm + "\t" + str(pos) + "\t")
					fout_sites_merged.write("\t".join([str(i) for i in m_sites_chrm[pos]]) + "\n")
				
				#this will use the number of clipped reads within the nearby region
				m_sites_chrm_filtered = xfilter.parse_sites_with_clip_cutoff_for_chrm(m_sites_chrm, cutoff_left_clip,
												      cutoff_right_clip,
												      cutoff_clip_mate_in_rep, cutoff_clip_mate_in_cns)
				for pos in m_sites_chrm_filtered:
					fout_sites.write(chrm + "\t" + str(pos) + "\t")
					fout_sites.write("\t".join([str(i) for i in m_sites_chrm_filtered[pos]]) + "\n")
				
				del m_sites_chrm_filtered # YW added to save memory
				del m_sites_chrm # YW added to save memory
				
		# YW 2021/04/07 added to remove loaded files
		for i in range(cnt):
			sf_tmp = self.working_folder + global_values.CLIP_TMP + "{0}".format(i)
			if os.path.isfile(sf_tmp):
				os.remove(sf_tmp)

####
	#this version is designed for mosaic calling only, which:
	#1. will count the polyA reads for each site
	#2. using polyA as a cutoff to filter out sites
	def call_TEI_candidate_sites_from_multiple_alignmts_mosaic(self, sf_annotation, sf_rep_cns, sf_ref, b_se,
								   cutoff_left_clip,
								   cutoff_right_clip, cutoff_clip_mate_in_rep, cutoff_polyA,
								   sf_clip_folder,
								   b_force, max_cov, sf_out):
		cnt = 0
		s_sample_bam = ""
		b_set = False
		with open(self.sf_list) as fin_bam_list:
			i_idx_bam = 0  # indicates which bam this is
			for line in fin_bam_list:  ###for each bam file
				fields = line.split()
				if len(fields) < 2:
					print("Error input bam file {0}!!!".format(line.rstrip()))
					continue
				sf_ori_bam = fields[0]
				s_read_type = fields[1].rstrip()
				print("Input bam {0} is sequenced from {1} platform!".format(sf_ori_bam, s_read_type))
				if global_values.X10 == s_read_type:  ###for 10X, we use a larger cutoff, as more clipped reads
					print("10X bam! Set the initial cutoff as {0}".format(global_values.INITIAL_MIN_CLIP_CUTOFF_X10))
					global_values.set_initial_min_clip_cutoff(global_values.INITIAL_MIN_CLIP_CUTOFF_X10)
				else:
					if cutoff_left_clip <= 3:#
						print("Clip cutoff is small (<=3), keep all the clipped reads (initial cutoff set as 1)!!!")
						global_values.set_initial_min_clip_cutoff(1)  # for low coverage data, set this to 1
					else:
						global_values.set_initial_min_clip_cutoff(global_values.INITIAL_MIN_CLIP_CUTOFF_ILLUMINA)
				if len(sf_ori_bam) <= 1:
					continue
				if b_set == False:
					s_sample_bam = sf_ori_bam
					b_set = True
				
				b_cutoff = True
				cutoff_hit_rep_copy = global_values.INITIAL_MIN_CLIP_CUTOFF
				
				# view the barcode bam as normal illumina bam
				# for each alignment, has one output
				sf_out_tmp = self.working_folder + global_values.CLIP_TMP + '{0}'.format(cnt)
				cnt += 1
				
				caller = TELocator(sf_ori_bam, sf_ori_bam, self.working_folder, self.n_jobs, self.sf_ref)
				# s_working_folder + global_values.CLIP_FOLDER + "/"+sf_bam_name + global_values.CLIP_FQ_SUFFIX
				sf_new_pub = ""
				if len(sf_clip_folder) == 0 or sf_clip_folder == None:
					print("public folder is null!!!!")
					continue
				if sf_clip_folder[-1] == "/":
					sf_new_pub = sf_clip_folder + "{0}/".format(i_idx_bam)
				else:
					sf_new_pub = sf_clip_folder + "/{0}/".format(i_idx_bam)
				caller.call_TEI_candidate_sites_from_clip_reads_v2_mosaic(sf_annotation, sf_rep_cns, sf_ref, b_se,
											  cutoff_hit_rep_copy, cutoff_hit_rep_copy,
											  b_cutoff, sf_new_pub, i_idx_bam, b_force, max_cov,
											  sf_out_tmp)
				i_idx_bam += 1
				####
		# get all the chromsomes names
		bam_info = BamInfo(s_sample_bam, self.sf_ref)
		b_with_chr = bam_info.is_chrm_contain_chr()
		m_chrms = bam_info.get_all_reference_names()
		
		xfilter = XIntermediateSites()
		xchrom = XChromosome()
		sf_out_merged = sf_out + "_tmp"
		with open(sf_out_merged, "w") as fout_sites_merged, open(sf_out, "w") as fout_sites:
			for chrm in m_chrms:  # write out chrm by chrm to save memory
				if xchrom.is_decoy_contig_chrms(chrm) == True:  ###filter out decoy and other contigs
					continue
				m_sites_chrm = {}
				for i in range(cnt):
					sf_tmp = self.working_folder + global_values.CLIP_TMP + "{0}".format(i)
					if os.path.isfile(sf_tmp) == False:
						print("Errors happen, file {0} doens't exist!".format(sf_tmp))
						continue
					with open(sf_tmp) as fin_tmp:
						for line in fin_tmp:
							fields = line.split()
							tmp_chrm = bam_info.process_chrm_name(fields[0], b_with_chr)
							if tmp_chrm != chrm:
								continue
							pos = int(fields[1])
							
							if pos not in m_sites_chrm:
								m_sites_chrm[pos] = []
								for value in fields[2:]:
									m_sites_chrm[pos].append(int(value))
							else:
								i_value = 0
								for value in fields[2:]:
									###sum (left-realign, right-realign, mate_in_rep, lpolyA, rpolyA)
									m_sites_chrm[pos][i_value] += int(value)
									i_value += 1
				
				for pos in m_sites_chrm:
					lth = len(m_sites_chrm[pos])
					fout_sites_merged.write(chrm + "\t" + str(pos) + "\t")
					for i in range(lth):
						s_feature = str(m_sites_chrm[pos][i])
						fout_sites_merged.write(s_feature + "\t")
					fout_sites_merged.write("\n")
				
				# this will use the number of clipped reads within the nearby region
				m_sites_chrm_filtered = xfilter.parse_sites_with_clip_cutoff_for_chrm_with_polyA(m_sites_chrm,
														 cutoff_left_clip,
														 cutoff_right_clip,
														 cutoff_clip_mate_in_rep,
														 cutoff_polyA)
				for pos in m_sites_chrm_filtered:
					lth = len(m_sites_chrm_filtered[pos])
					fout_sites.write(chrm + "\t" + str(pos) + "\t")
					for i in range(lth):
						s_feature = str(m_sites_chrm_filtered[pos][i])
						fout_sites.write(s_feature + "\t")
					fout_sites.write("\n")

####
	# For given candidate sites from clip reads,
	# sum the num of the discordant pairs from different alignments
	# YW 2020/08/03 github update, 2 more arguments
	# YW 2021/04/23 change r_lcluster, s_lc_chrm, s_lc_pos, r_rcluster, s_rc_chrm, s_rc_pos from list to a single entry
	# YW 2021/04/27 add f_disc_fa_tmp for extract_mate_reads_by_name function
	# YW 2021/04/30 remove sf_out from input argument
	# YW 2021/05/19 removed b_tumor
	def filter_candidate_sites_by_discordant_pairs_multi_alignmts(self, m_sites, iext, i_is, f_dev, cutoff,
									sf_rep_cns_Alu, sf_rep_cns_L1, sf_rep_cns_SVA,
									sf_annotation_Alu, sf_annotation_L1, sf_annotation_SVA, 
									sf_raw_disc=""):
		sf_disc_working_folder = self.working_folder + global_values.DISC_FOLDER
		cmd = f"mkdir -p {sf_disc_working_folder}"
		self.cmd_runner.run_cmd_small_output(cmd)
		sf_disc_working_folder += '/'
		m_disc = {}
		with open(self.sf_list) as fin_list:
			cnt = 0
			for line in fin_list:
				fields=line.split()
				sf_bam = fields[0]
                
				caller = TELocator(sf_bam, sf_bam, self.working_folder, self.n_jobs, self.sf_ref)
				# YW 2020/08/01 added the global_value (originally just 1)
				tmp_cutoff = global_values.INITIAL_MIN_DISC_CUTOFF
				#############################################################
				# YW 2020/07/04 added the following if statement
				if cutoff <= 1:
					tmp_cutoff = 0
					print("Disc cutoff is small (<=1) , we are using {0} for initial cutoff".format(tmp_cutoff))
					global_values.set_initial_min_disc_cutoff(tmp_cutoff)
				
				else:
					print("We are using {0} for initial cutoff".format(tmp_cutoff))
				####################################################
				# 2021/04/27 YW makes temporary directory to store disc_pos by sample
				tmp_pos_folder = sf_disc_working_folder + str(cnt) + "/"
				cmd = f"mkdir -p {tmp_pos_folder}"
				self.cmd_runner.run_cmd_small_output(cmd)
				# YW 2021/04/27 add disc_pos output
				sf_disc_pos_tmp = "".join([tmp_pos_folder, "sample", str(cnt), global_values.DISC_POS_SUFFIX]) # YW 2021/04/26 moved from cns-remapping --> change this
				####################################################
				# 2020/08/03 github update: more output for filter_candidate_sites_by_discordant_pairs_non_barcode
				# m_sites_discord, m_sites_raw_disc = caller.filter_candidate_sites_by_discordant_pairs_non_barcode(
				# m_sites, iext, i_is, f_dev, sf_annotation_Alu, sf_annotation_L1, sf_annotation_SVA, tmp_cutoff, sf_disc_pos_tmp, sf_disc_working_folder, tmp_pos_folder) # YW 2021/04/27 added sf_disc_pos_tmp, sf_disc_working_folder, tmp_pos_folder variables, YW TO DO: remove m_sites_discord from output
				m_sites_raw_disc = caller.filter_candidate_sites_by_discordant_pairs_non_barcode_short(
							m_sites, iext, i_is, f_dev, sf_annotation_Alu, sf_annotation_L1, sf_annotation_SVA, tmp_cutoff, sf_disc_pos_tmp, sf_disc_working_folder, tmp_pos_folder) # YW 2021/04/27 added sf_disc_pos_tmp, sf_disc_working_folder, tmp_pos_folder variables, YW 2021/04/29 removed m_sites_discord from output
				##########################################################
				# YW 2021/04/27 added to extract disc reads and realign to repeat cns
				bam_info = BamInfo(sf_bam, self.sf_ref)
				sf_disc_fa_tmp = tmp_pos_folder + "temp_disc.fa"
				bam_info.extract_mate_reads_by_name_short(sf_disc_pos_tmp, global_values.BIN_SIZE, sf_disc_working_folder, self.n_jobs, sf_disc_fa_tmp)
				##########################################################
				
				##########################################################
				# YW 2021/04/27 copied and modified from from x_clip_disc_filter.py to realign and parse disc reads (now doing these by sample, not all at once)
				# YW 2021/05/26 made this chunk parallel
				sf_disc_algnmt_Alu = tmp_pos_folder + "Alu" + global_values.DISC_SAM_SUFFIX
				sf_disc_algnmt_L1 = tmp_pos_folder + "L1" + global_values.DISC_SAM_SUFFIX
				sf_disc_algnmt_SVA = tmp_pos_folder + "SVA" + global_values.DISC_SAM_SUFFIX
				bwa_align = BWAlign(global_values.BWA_PATH, global_values.BWA_REALIGN_CUTOFF, self.n_jobs)
				
				# bwa_align.realign_disc_reads(sf_rep_cns_Alu, sf_disc_fa_tmp, sf_disc_algnmt_Alu)
				# bwa_align.realign_disc_reads(sf_rep_cns_L1, sf_disc_fa_tmp, sf_disc_algnmt_L1)
				# bwa_align.realign_disc_reads(sf_rep_cns_SVA, sf_disc_fa_tmp, sf_disc_algnmt_SVA)
				
				#######################################
				# to parallel disc realignment
				cns_parallel = CNS_PARALLEL(self.n_jobs, self.working_folder, sf_disc_fa_tmp, cnt)
				L1_script, L1_done, L1_fail = cns_parallel.gnrt_disc_py_scripts("L1", sf_rep_cns_L1, sf_disc_algnmt_L1)
				SVA_script, SVA_done, SVA_fail = cns_parallel.gnrt_disc_py_scripts("SVA", sf_rep_cns_SVA, sf_disc_algnmt_SVA)
				cns_parallel.run_sbatch_scripts([(L1_script, L1_done, L1_fail), (SVA_script, SVA_done, SVA_fail)], "disc")
				
				bwa_align.realign_disc_reads(sf_rep_cns_Alu, sf_disc_fa_tmp, sf_disc_algnmt_Alu)
				print("Alu cns realignment has finished!")
				while (not os.path.exists(L1_done) or not os.path.exists(SVA_done)) and (not os.path.exists(L1_fail) or not os.path.exists(SVA_fail)):
					sleep(global_values.CHECK_INTERVAL)
				if os.path.exists(L1_fail):
					sys.exit("L1 cns realignment failed.")
				if os.path.exists(SVA_fail):
					sys.exit("SVA cns realignment failed.")
				# YW 2021/05/26 remove files so next time we run this chunk of code we don't automatically skip L1 and SVA realignment
				os.remove(L1_done)
				os.remove(SVA_done)
				if os.path.exists(sf_disc_algnmt_L1) and os.path.exists(sf_disc_algnmt_SVA):
					if os.path.getsize(sf_disc_algnmt_L1)>0 and os.path.getsize(sf_disc_algnmt_SVA)>0:
						print(f"Alu, L1, SVA disc realignments have all finished for bam {cnt}. Proceed to counting clipped parts...")
					else:
						sys.exit("Something went wrong with L1 or SVA disc realignment!!!")
				else:
					sys.exit("Something went wrong with L1 or SVA disc realignment!!!")
				########################################
				
				os.remove(sf_disc_fa_tmp)
				self.parse_disc_algnmt_consensus_short(sf_disc_algnmt_Alu, global_values.BMAPPED_CUTOFF, m_disc, "Alu")
				os.remove(sf_disc_algnmt_Alu)
				self.parse_disc_algnmt_consensus_short(sf_disc_algnmt_L1, global_values.BMAPPED_CUTOFF, m_disc, "L1")
				os.remove(sf_disc_algnmt_L1)
				self.parse_disc_algnmt_consensus_short(sf_disc_algnmt_SVA, global_values.BMAPPED_CUTOFF, m_disc, "SVA")
				os.remove(sf_disc_algnmt_SVA)
				##########################################################
				
				xfilter = XIntermediateSites()
				# sf_out_tmp = self.working_folder + global_values.DISC_TMP + str(cnt)
				# xfilter.output_candidate_sites(m_sites_discord, sf_out_tmp) # YW 2021/04/28 commented out
				# 2020/08/03 github update: add the following 2 lines
				sf_raw_tmp=self.working_folder+global_values.RAW_DISC_TMP+str(cnt)
				xfilter.output_candidate_sites(m_sites_raw_disc, sf_raw_tmp)
				cnt += 1
		 
####		# YW 2020/08/03 github update (all the following code)
		if sf_raw_disc=="":
			return
		#####
		# merge the output by summing up all the alignments,
		#  and output in a single file
		m_merged_raw_sites = {}
		for i in range(cnt):
			sf_tmp = self.working_folder + global_values.RAW_DISC_TMP + str(i)
			with open(sf_tmp) as fin_tmp:
				for line in fin_tmp:
					fields = line.split()
					chrm = fields[0]
					pos = int(fields[1])
					n_raw_left_disc = int(fields[2])
					n_raw_right_disc = int(fields[3])
					n_left_disc_Alu = int(fields[4])
					n_right_disc_Alu = int(fields[5])
					n_left_disc_L1 = int(fields[6])
					n_right_disc_L1 = int(fields[7])
					n_left_disc_SVA = int(fields[8])
					n_right_disc_SVA = int(fields[9])
					r_lcluster=float(fields[10])
					# s_lc_chrm=fields[11]
					# s_lc_pos=fields[12]
					r_rcluster=float(fields[13])
					# s_rc_chrm=fields[14]
					# s_rc_pos=fields[15]
					
					if chrm not in m_merged_raw_sites:
						m_merged_raw_sites[chrm] = {}
					if pos not in m_merged_raw_sites[chrm]:
						m_merged_raw_sites[chrm][pos] = [n_raw_left_disc, n_raw_right_disc,
										 n_left_disc_Alu, n_right_disc_Alu,
										 n_left_disc_L1, n_right_disc_L1,
										 n_left_disc_SVA, n_right_disc_SVA,
										 r_lcluster, r_rcluster]
					else:
						# YW added below because all samples only have one bam file
						sys.exit(f"{chrm} {pos} already exists in dict m_merged_raw_sites")
						# m_merged_raw_sites[chrm][pos][0] += n_raw_left_disc
						# m_merged_raw_sites[chrm][pos][1] += n_raw_right_disc
						# m_merged_raw_sites[chrm][pos][2] += n_left_disc_Alu
						# m_merged_raw_sites[chrm][pos][3] += n_right_disc_Alu
						# m_merged_raw_sites[chrm][pos][4] += n_left_disc_L1
						# m_merged_raw_sites[chrm][pos][5] += n_right_disc_L1
						# m_merged_raw_sites[chrm][pos][6] += n_left_disc_SVA
						# m_merged_raw_sites[chrm][pos][7] += n_right_disc_SVA
						# m_merged_raw_sites[chrm][pos][8].append(r_lcluster) # YW: may wanna get max
						# m_merged_raw_sites[chrm][pos][9].append(s_lc_chrm) # YW: may wanna see if they all match to the same chr
						# m_merged_raw_sites[chrm][pos][10].append(s_lc_pos) # YW: may wanna get most frequent
						# m_merged_raw_sites[chrm][pos][11].append(r_rcluster) # YW: may wanna get max
						# m_merged_raw_sites[chrm][pos][12].append(s_rc_chrm) # YW: may wanna see if they all match to the same chr
						# m_merged_raw_sites[chrm][pos][13].append(s_rc_pos) # YW: may wanna get most frequent
						#####################################################
						# YW 2021/04/29 added below to output counts of disc reads aligned to TE cns
			os.remove(sf_tmp) # YW 2021/04/07 added to remove merged file
			
		# YW 2021/05/10 add disc reads mapping to cns counts
		for chrm in m_merged_raw_sites:
			if chrm not in m_disc:
				for pos in m_merged_raw_sites[chrm]:
					m_merged_raw_sites[chrm][pos].extend([0, 0, 0])
			else:
				for pos in m_merged_raw_sites[chrm]:
					if pos in m_disc[chrm]:
						m_merged_raw_sites[chrm][pos].extend(m_disc[chrm][pos])
					else:
						m_merged_raw_sites[chrm][pos].extend([0, 0, 0])
		#####################################################
                
		i_half_cutoff=cutoff/2
		# if b_tumor==False:#for germline, set a little bit higher cutoff
		# 	i_half_cutoff= int(cutoff*3/4) + 1
		if i_half_cutoff<1:
			# YW 2020/08/11 changed this
			# i_half_cutoff=1
			i_half_cutoff=0
		
		# with open(sf_raw_disc, "w") as fout_sites:
		#     for chrm in m_merged_raw_sites:
		#         for pos in m_merged_raw_sites[chrm]:
		#             n_raw_left = m_merged_raw_sites[chrm][pos][0]
		#             n_raw_right = m_merged_raw_sites[chrm][pos][1]
		#             n_left = m_merged_raw_sites[chrm][pos][2]#fall in repetitive region
		#             n_right = m_merged_raw_sites[chrm][pos][3]#fall in repetitive region
		# 
		#             #left consistent
		#             b_l_consistent=False
		#             #here "1" indicates True, which means
		#             if (n_raw_left>=i_half_cutoff) and ("1" in m_merged_raw_sites[chrm][pos][4]):
		#                 b_l_consistent=True
		#             #right consistent
		#             b_r_consistent=False
		#             if (n_raw_right >= i_half_cutoff) and ("1" in m_merged_raw_sites[chrm][pos][7]):
		#                 b_r_consistent=True
		#             if b_l_consistent or b_r_consistent:
		#                 fout_sites.write(
		#                     chrm + "\t" + str(pos) + "\t" + str(n_raw_left) + "\t" + str(n_raw_right) + "\t"
		#                     + str(n_left) + "\t" + str(n_right) + "\n")
		# YW 2021/04/17 changed the lines above to output all the information in m_merged_raw_sites
		# YW 2021/04/22 fixed bugs
		with open(sf_raw_disc, "w") as fout_sites:
			for chrm in m_merged_raw_sites:
				for pos in m_merged_raw_sites[chrm]:
					fout_sites.write("\t".join([chrm, str(pos), ""]))
					new_lst = list(map(str, m_merged_raw_sites[chrm][pos]))
					fout_sites.write("\t".join(new_lst) + "\n")
####
####
	# This function is not used here
	# merge two dictionary
	def merge_sites(self, m_tmp, m_final):
		for chrm in m_tmp:
			if chrm not in m_final:
				m_final[chrm] = {}
			for pos in m_tmp[chrm]:
				if pos not in m_final[chrm]:
					m_final[chrm][pos] = m_tmp[chrm][pos]
				else:
					lth = len(m_tmp[chrm][pos])
					for i in range(lth):
						m_final[chrm][pos][i] += m_tmp[chrm][pos][i]
	
	###################################################################
	# YW 2021/04/28 copied from x_clip_disc_filter.py
	# used to be parse_disc_algnmt_consensus
	# remove sample_id from m_disc_sample dict (now doing this for every single sample)
	# get rid of irrelevant info
	# add m_disc_sample as input, rather than returning a new dict
	# add rep_type
	# TO DO: RETURN THE MAPPED DISC READS
	def parse_disc_algnmt_consensus_short(self, sf_disc_alignmt, bmapped_cutoff, m_disc_sample, rep_type):
		samfile = pysam.AlignmentFile(sf_disc_alignmt, "r", reference_filename=self.sf_ref)
		for algnmt in samfile.fetch():
			####also, for clipped mapped reads, need to check whether the clipped parts can be split to two parts!!!!!!
			# YW 2020/07/13 updated the format comment
			# fmt: read_id~b_lclip~b_rclip~is_anchor_rc~is_anchor_mate_rc~anchor_map_pos~s_insertion_chrm~s_insertion
			# YW 2021/04/28 don't need this much info, change fmt to read_id~anchor_map_pos~s_insertion_chrm~s_insertion
			# YW 2020/08/03 github update: added the following chunk
			if algnmt.is_unmapped == True:  ####skip the unmapped reads
				continue
			# first check whether read is qualified mapped
			l_cigar = algnmt.cigar
			read_seq = algnmt.query_sequence # YW 2020/08/15 moved up here ==> return read_seq later
			# YW 2020/08/15 changed the input of is_clipped_part_qualified_algnmt to consider mask
			b_clip_qualified_algned = self.is_clipped_part_qualified_algnmt_disc(l_cigar, bmapped_cutoff, read_seq)
			if b_clip_qualified_algned == False:  # skip the unqualified re-aligned parts
				continue

			read_info = algnmt.query_name
			# read_seq = algnmt.query_sequence # YW TO DO!
			read_info_fields = read_info.split(global_values.SEPARATOR)
			# reads info = >read_name~is_first~anchor_map_pos~insertion_chrm~insertion_pos
			# s_anchor_lclip=read_info_fields[-7] #anchor read is left clip or not: 1 indicates clip
			# s_anchor_rclip=read_info_fields[-6] #anchor read is right clip or not: 1 indicates clip
			# s_anchor_rc = read_info_fields[-5]
			# s_anchor_mate_rc = read_info_fields[-4]
			anchor_map_pos = int(read_info_fields[-3]) #this is the position of the left-most mapped base
			ins_chrm = read_info_fields[-2]
			ins_pos = int(read_info_fields[-1])

			if ins_chrm not in m_disc_sample:
				m_disc_sample[ins_chrm] = {}
			if ins_pos not in m_disc_sample[ins_chrm]:
				m_disc_sample[ins_chrm][ins_pos] = [0, 0, 0]
			if rep_type == "Alu":
				m_disc_sample[ins_chrm][ins_pos][0] += 1
			elif rep_type == "L1":
				m_disc_sample[ins_chrm][ins_pos][1] += 1
			elif rep_type == "SVA":
				m_disc_sample[ins_chrm][ins_pos][2] += 1
		
		samfile.close()
		# return m_disc_sample
		
	# check the clipped part is qualified aligned or not
	# YW 2020/08/15 changed the input of is_clipped_part_qualified_algnmt to consider mask 2 (not counting them in the clip_len for clipped and n_total)
	# used to be is_clipped_part_qualified_algnmt
	# 2021/05/03 YW got rid of n_map output 
	def is_clipped_part_qualified_algnmt_disc(self, l_cigar, ratio_cutoff, seq): # YW added seq argument
		if len(l_cigar) < 1:  # wrong alignment
			#return False, 0
			return False
		if len(l_cigar) > 2:
			####check the cigar
			###if both clipped, and the clipped part is large, then skip
			b_left_clip = False
			i_left_clip_len = 0
			if l_cigar[0][0] == 4 or l_cigar[0][0] == 5:  # left clipped
				b_left_clip = True
				i_left_clip_len = l_cigar[0][1]
				# YW 2020/08/15 added the following 4 lines
				i_left_clip_seq = seq[:l_cigar[0][1]]
				for ch in i_left_clip_seq:
					if ch == "N":
						i_left_clip_len -= 1
			b_right_clip = False
			i_right_clip_len = 0
			if l_cigar[-1][0] == 4 or l_cigar[-1][0] == 5:  # right clipped
				b_right_clip = True
				i_right_clip_len = l_cigar[-1][1]
				# YW 2020/08/15 added the following 5 lines
				start_pos = -1 * l_cigar[-1][1]
				i_right_clip_seq = seq[start_pos:]
				for ch in i_right_clip_seq:
					if ch == "N":
						i_right_clip_len -= 1

			if b_left_clip == True and b_right_clip == True:
				if (i_left_clip_len > global_values.MAX_CLIP_CLIP_LEN) and (i_right_clip_len > global_values.MAX_CLIP_CLIP_LEN):
					#return False, 0
					return False
		
		####YW updated the comment 2020/07/11: for the alignment (of the clipped read), if the mapped part is smaller than bmapped_cutoff * clipped part,
		####then skip
		n_total = 0
		n_map = 0
		for (type, lenth) in l_cigar:
			if type == 0:
				n_map += lenth
			if type != 2:  # deletion is not added to the total length
				n_total += lenth
		# YW 2020/08/15 added the following 3 lines
		for ch in seq:
			if ch == "N":
				n_total -= 1
		if n_map < (n_total * ratio_cutoff):  ########################require at least bmapped_cutoff of the seq is mapped !!!!!!!!
			#return False, 0
			return False
		#return True, n_map
		return True


def unwrap_self_filter_by_discordant(arg, **kwarg):
	return TELocator.run_filter_by_discordant_pair_by_chrom(*arg, **kwarg)


def unwrap_self_filter_by_discordant_non_barcode(arg, **kwarg):
	return TELocator.run_filter_by_discordant_pair_by_chrom_non_barcode(*arg, **kwarg)

def unwrap_self_filter_by_discordant_non_barcode_short(arg, **kwarg):
	return TELocator.run_filter_by_discordant_pair_by_chrom_non_barcode_short(*arg, **kwarg)


def unwrap_self_filter_by_barcode_coverage(arg, **kwarg):
	return TELocator.run_filter_by_barcode_coverage(*arg, **kwarg)


class TELocator():
	def __init__(self, sf_bam, sf_barcode_bam, s_working_folder, n_jobs, sf_ref):
		self.sf_bam = sf_bam
		self.sf_barcode_bam = sf_barcode_bam
		self.working_folder = s_working_folder
		self.n_jobs = int(n_jobs)
		self.sf_reference = sf_ref  ##reference genome
		self.cmd_runner = CMD_RUNNER()

	###First, Use (left, right) clipped read as threshold. Also, require some of the mate read are within repeat region
	###Then: check the nearby small region, whether the merged number satisfy the threshold
	###Then, from the candidate list, pick the peak in each window.
	
	###First, Use (left, right) clipped read as threshold. Also, require some of the mate read are within repeat region
	##Note, this version consider the insertion with deletion cases, that is common in many cases
 	##So, for TEI with deletion, there will be two breakpoints, and at each breakpoint, only one type of clipped reads
 	# YW 2020/08/09: this function doesn't require mate read in repeat region (--cr), the filtering is done in the last step after merging nearby clipping positions later on
	# YW 2020/08/09: did not see how deletion is considered.
	def call_TEI_candidate_sites_from_clip_reads_v2(self, sf_annotation_Alu, sf_annotation_L1, sf_annotation_SVA,
							sf_rep_cns_Alu, sf_rep_cns_L1, sf_rep_cns_SVA,
							sf_rep_Alu, sf_rep_L1, sf_rep_SVA,
							b_se, cutoff_left_clip,
							cutoff_right_clip, b_cutoff, sf_pub_folder, idx_bam,
							b_force, max_cov_cutoff, sf_out): # YW 2021/03/18 take out max_cov_cutoff later
		# this is a public folder for different type of repeats to share the clipped reads
		if sf_pub_folder[-1]!="/":
			sf_pub_folder+="/"
		if os.path.exists(sf_pub_folder) == False:
			cmd = "mkdir -p {0}".format(sf_pub_folder)
			self.cmd_runner.run_cmd_small_output(cmd)
		
		#this is the local folder for the current read type to save the tmp files
		sf_clip_working_folder = self.working_folder + global_values.CLIP_FOLDER + "/{0}/".format(idx_bam)
		if os.path.exists(sf_clip_working_folder) == False:
			cmd = "mkdir -p {0}".format(sf_clip_working_folder)
			self.cmd_runner.run_cmd_small_output(cmd)
		
		clip_info = ClipReadInfo(self.sf_bam, self.n_jobs, self.sf_reference)
		######1. so first, re-align the clipped parts, and count the number of supported clipped reads
		####gnrt the clipped parts file
		sf_bam_name = os.path.basename(self.sf_bam)
		sf_all_clip_fq = sf_pub_folder + sf_bam_name + global_values.CLIP_FQ_SUFFIX
		clip_info.set_working_folder(sf_clip_working_folder)
		if os.path.islink(sf_all_clip_fq)==False or b_force==True:
			print("Collected clipped reads file {0} doesn't exist. Generate it now!".format(sf_all_clip_fq))
			##collect the clip positions
			initial_clip_pos_freq_cutoff = global_values.INITIAL_MIN_CLIP_CUTOFF ##########################################################################
			print("Initial minimum clip cutoff is {0}".format(initial_clip_pos_freq_cutoff))
			clip_info.collect_clip_positions(sf_annotation_Alu, sf_annotation_L1, sf_annotation_SVA,
							initial_clip_pos_freq_cutoff, b_se, sf_pub_folder) ##save clip pos by chrm
			print("Output info: Collect clipped parts for file ", self.sf_bam)
			sf_all_clip_fq_ori=sf_clip_working_folder+sf_bam_name + global_values.CLIP_FQ_SUFFIX
			clip_info.collect_clipped_parts(sf_all_clip_fq_ori)
			
			if os.path.isfile(sf_all_clip_fq)==True or os.path.islink(sf_all_clip_fq)==True:
				os.remove(sf_all_clip_fq)
			cmd="ln -s {0} {1}".format(sf_all_clip_fq_ori, sf_all_clip_fq)
			self.cmd_runner.run_cmd_small_output(cmd)
		else:
			print("Collected clipped reads file {0} already exist!".format(sf_all_clip_fq))
####
		####align the clipped parts to repeat copies
		# YW 2021/03/18 add Alu, L1, SVA
		sf_algnmt_Alu = self.working_folder + sf_bam_name + global_values.CLIP_BAM_SUFFIX + ".Alu"
		sf_algnmt_L1 = self.working_folder + sf_bam_name + global_values.CLIP_BAM_SUFFIX + ".L1"
		sf_algnmt_SVA = self.working_folder + sf_bam_name + global_values.CLIP_BAM_SUFFIX + ".SVA"
		print("Output info: Re-align clipped parts for file ", self.sf_bam)
		
		# YW 2021/05/26 modified the following to parallelize clipped read mapping to cns
		bwa_align = BWAlign(global_values.BWA_PATH, global_values.BWA_REALIGN_CUTOFF, self.n_jobs)
		cns_parallel = CNS_PARALLEL(self.n_jobs, self.working_folder, sf_all_clip_fq, idx_bam)
		L1_script, L1_done, L1_fail = cns_parallel.gnrt_clip_py_scripts("L1", sf_rep_cns_L1, sf_rep_L1, sf_algnmt_L1)
		SVA_script, SVA_done, SVA_fail = cns_parallel.gnrt_clip_py_scripts("SVA", sf_rep_cns_SVA, sf_rep_SVA, sf_algnmt_SVA)
		cns_parallel.run_sbatch_scripts([(L1_script, L1_done, L1_fail), (SVA_script, SVA_done, SVA_fail)], "clip")
		# YW 2021/03/18 add Alu, L1, SVA
		bwa_align.two_stage_realign(sf_rep_cns_Alu, sf_rep_Alu, sf_all_clip_fq, sf_algnmt_Alu)
		# bwa_align.two_stage_realign(sf_rep_cns_L1, sf_rep_L1, sf_all_clip_fq, sf_algnmt_L1)
		# bwa_align.two_stage_realign(sf_rep_cns_SVA, sf_rep_SVA, sf_all_clip_fq, sf_algnmt_SVA)
		print("Alu cns realignment has finished!")
		while not os.path.exists(L1_done) or not os.path.exists(SVA_done) and (not os.path.exists(L1_fail) or not os.path.exists(SVA_fail)):
			sleep(global_values.CHECK_INTERVAL)
		if os.path.exists(L1_fail):
			sys.exit("L1 cns realignment failed.")
		if os.path.exists(SVA_fail):
			sys.exit("SVA cns realignment failed.")
		# YW 2021/05/26 remove files so next time we run this chunk of code we don't automatically skip L1 and SVA realignment
		os.remove(L1_done)
		os.remove(SVA_done)
		if os.path.exists(sf_algnmt_L1) and os.path.exists(sf_algnmt_SVA):
			if os.path.getsize(sf_algnmt_L1)>0 and os.path.getsize(sf_algnmt_SVA)>0:
				print(f"Alu, L1, SVA clip realignments have all finished for bam {idx_bam}. Proceed to counting clipped parts...")
			else:
				sys.exit("Something went wrong with L1 or SVA clip realignment!!!")
		else:
			sys.exit("Something went wrong with L1 or SVA clip realignment!!!")
		# remove unused intermediate files to save disk space
		os.remove(sf_all_clip_fq)
		
		####cnt number of clipped reads aligned to repeat copies from the re-alignment
		# YW 2021/03/18 add Alu, L1, SVA
		clip_info.cnt_clip_part_aligned_to_rep(sf_algnmt_Alu, sf_algnmt_L1, sf_algnmt_SVA)  ##require at least 3/4 (YW changed from half from previous comment) of the seq is mapped !!!!

		# if b_cutoff is set, then directly return the dict
		if b_cutoff == False:
			clip_info.merge_clip_positions(sf_pub_folder, sf_out)
		else:
			clip_info.merge_clip_positions_with_cutoff(cutoff_left_clip, cutoff_right_clip, max_cov_cutoff, sf_pub_folder, sf_out)
		if os.path.isfile(sf_algnmt_Alu)==True:####remove the file
			os.remove(sf_algnmt_Alu)
		if os.path.isfile(sf_algnmt_L1)==True:####remove the file
			os.remove(sf_algnmt_L1)
		if os.path.isfile(sf_algnmt_SVA)==True:####remove the file
			os.remove(sf_algnmt_SVA)
####

	def collect_all_clipped_reads_only(self, sf_annotation, b_se, s_working_folder, wfolder_pub_clip):
		sf_clip_working_folder = s_working_folder + global_values.CLIP_FOLDER + "/"
		if len(sf_clip_working_folder)>1 and sf_clip_working_folder[-1]!="/":
			sf_clip_working_folder+="/"
		if os.path.exists(sf_clip_working_folder) == False:
			cmd = "mkdir {0}".format(sf_clip_working_folder)
			self.cmd_runner.run_cmd_small_output(cmd)
		
		clip_info = ClipReadInfo(self.sf_bam, self.n_jobs, self.sf_reference)
		clip_info.set_working_folder(sf_clip_working_folder)
		
		##collect the clip positions
		##in format {chrm: {map_pos: [left_cnt, right_cnt, mate_within_rep_cnt]}}
		print("Output info: Collect clip positions for file ", self.sf_bam)
		initial_clip_pos_freq_cutoff = 2  ##############################################################################
		##save clip pos by chrm
		clip_info.collect_clip_positions(sf_annotation, initial_clip_pos_freq_cutoff, b_se, wfolder_pub_clip)
		
		######1. so first, re-align the clipped parts, and count the number of supported clipped reads
		####gnrt the clipped parts file
		sf_bam_name = os.path.basename(self.sf_bam)
		sf_all_clip_fq = sf_clip_working_folder + sf_bam_name + global_values.CLIP_FQ_SUFFIX
		print("Output info: Collect clipped parts for file ", self.sf_bam)
		# if os.path.isfile(sf_all_clip_fq)==False:
		clip_info.collect_clipped_parts(sf_all_clip_fq)
####

	def _is_decoy_contig_chrms(self, chrm):
		fields = chrm.split("_")
		if len(fields) > 1:
			return True
		elif chrm == "hs37d5":
			return True
		else:
			return False
		
	###Input in format {chrm: {map_pos: [left_cnt, right_cnt, mate_within_rep_cnt]}}
	def filter_out_decoy_contig_chrms(self, m_candidate_list):
		m_new_list = {}
		for chrm in m_candidate_list:
			if self._is_decoy_contig_chrms(chrm) == True:
				continue
			
			if chrm not in m_new_list:
				m_new_list[chrm] = {}
			for pos in m_candidate_list[chrm]:
				if pos not in m_new_list[chrm]:
					m_new_list[chrm][pos] = []
				for value in m_candidate_list[chrm][pos]:
					m_new_list[chrm][pos].append(value)
		return m_new_list


	#####First, Use (left, right) clipped read as threshold. Also, require some of the mate read are within repeat region
	##Note: this version will miss some cases, like insertion with deletion ones
	def call_TEI_candidate_sites_from_clip_reads(self, sf_ref, sf_annotation, cutoff_left_clip, cutoff_right_clip):
		clip_info = ClipReadInfo(self.sf_bam, self.n_jobs, self.sf_reference)
		clip_info.set_working_folder(self.working_folder)
		sf_all_clip_fq = self.working_folder + global_values.CLIP_FQ_SUFFIX
		m_clip_pos_freq = clip_info.collect_clipped_reads_with_position(sf_all_clip_fq, sf_annotation)
		
		# xself.output_candidate_sites(m_clip_pos_freq, "initial_clip_pos.txt")
		####Here need to use the cutoff to remove most of the unnecessary sites first !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		clip_info.correct_num_of_clip_by_realignment(sf_ref, sf_annotation, sf_all_clip_fq, m_clip_pos_freq)
		mate_in_rep_cutoff = (cutoff_right_clip + cutoff_left_clip) / 2
		m_candidate_sites = {}
		for chrm in m_clip_pos_freq:
			for pos in m_clip_pos_freq[chrm]:
				####here need to check the nearby region
				nearby_left_freq = 0
				nearby_right_freq = 0
				nearby_mate_in_rep = 0
				for i in range(-1 * global_values.NEARBY_REGION, global_values.NEARBY_REGION):
					i_tmp_pos = pos + i
					if i_tmp_pos in m_clip_pos_freq[chrm]:
						nearby_left_freq += m_clip_pos_freq[chrm][i_tmp_pos][0]
						nearby_right_freq += m_clip_pos_freq[chrm][i_tmp_pos][1]
						nearby_mate_in_rep += m_clip_pos_freq[chrm][i_tmp_pos][2]
				
				if nearby_left_freq >= cutoff_left_clip and nearby_right_freq >= cutoff_right_clip \
					and nearby_mate_in_rep >= mate_in_rep_cutoff:
					if chrm not in m_candidate_sites:
						m_candidate_sites[chrm] = {}
						# if pos not in m_candidate_sites[chrm]:
						#    m_candidate_sites[chrm] = {}
						i_left_cnt = m_clip_pos_freq[chrm][pos][0]
						i_right_cnt = m_clip_pos_freq[chrm][pos][1]
						i_mate_in_rep_cnt = m_clip_pos_freq[chrm][pos][2]
						m_candidate_sites[chrm][pos] = (i_left_cnt, i_right_cnt, i_mate_in_rep_cnt)
		return m_candidate_sites

    # #align the reads to the repeat copies to collect the repeat related reads only
    # #note, we also include two flank regions of the repeats to help find the transductions
    # def collect_reads_fall_in_repeats(self, m_site_algnmts, m_selected_algnmts, sf_out):
    #     #first, dump all the alignmnts to a file
    #     #then,

	def run_filter_by_discordant_pair_by_chrom(self, record):
		site_chrm1 = record[0]
		sf_bam = record[1]
		sf_barcode_bam = record[2]
		iextend = int(record[3])  ###extend some region on both sides in order to collect all barcodes
		i_is = int(record[4])
		f_dev = int(record[5])
		sf_annotation = record[6]
		sf_disc_working_folder = record[7]
		s_suffix = record[8]
		# iextend_small=record[9] ##extend a small region to compare the barcode difference
		# iextend_small = 300  ##extend a small region to compare the barcode difference ###############################
		
		sf_candidate_list = sf_disc_working_folder + site_chrm1 + s_suffix
		if os.path.exists(sf_candidate_list) == False:
			return
		m_candidate_pos = {}
		with open(sf_candidate_list) as fin_list:
			for line in fin_list:
				fields = line.split()
				pos = int(fields[1])
				m_candidate_pos[pos] = "\t".join(fields[2:])
				
		bam_info = BamInfo(sf_bam, self.sf_reference)
		b_with_chr = bam_info.is_chrm_contain_chr()
		m_chrms = bam_info.get_all_reference_names()
		site_chrm = bam_info.process_chrm_name(site_chrm1, b_with_chr)
		if site_chrm not in m_chrms:
			return
		xannotation = XAnnotation(sf_annotation)
		xannotation.set_with_chr(b_with_chr)
		xannotation.load_rmsk_annotation()
		xannotation.index_rmsk_annotation()
		
		m_new_candidate_sites = {}
		xbam = XBamInfo(sf_bam, sf_barcode_bam, self.sf_reference)
		xbam.index_reference_name_id()
		bamfile = xbam.open_bam_file(sf_bam)  ##open bam file
		barcode_bamfile = xbam.open_bam_file(sf_barcode_bam)  ##open barcode bam file
		for site_pos in m_candidate_pos:  ####candidate site position # structure: {barcode:[alignmts]}
			if site_pos < iextend:
				continue
			# n_barcode_diff, n_barcode_share = xbam.check_barcode_diff_v2(bamfile, site_chrm, site_pos, iextend)
			# print(site_chrm1, site_pos, " test1!!!!!!!!!")#############################################################
			m_site_algnmts = xbam.parse_alignments_for_one_site_v2(bamfile, barcode_bamfile, site_chrm, site_pos, iextend)
			if m_site_algnmts == None:  ####!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				continue
			n_barcode = len(m_site_algnmts)
			n_discordant_pair = 0
			n_one_in_rep_region_pair = 0
			n_both_in_rep_region_pair = 0
			for barcode in m_site_algnmts:
				l_algnmts = m_site_algnmts[barcode]
				for algnmt in l_algnmts:  ####Note, alignment here is in converted format: barcode as chromosome
					barcode_algnmt = BarcodeAlignment(algnmt)
					chrm = barcode_algnmt.get_chrm()
					map_pos = barcode_algnmt.get_map_pos()
					mate_chrm = barcode_algnmt.get_mate_chrm()
					mate_map_pos = int(algnmt.next_reference_start)
					template_lth = algnmt.template_length
					
					xalgnmt = XAlignment()
					b_discordant = xalgnmt.is_discordant_pair_no_unmap(chrm, mate_chrm, template_lth, i_is, f_dev)
					if b_discordant == True:
						n_discordant_pair += 1
					if algnmt.is_unmapped or algnmt.mate_is_unmapped:
						continue
					
					###1. consider first/mate reads only once?
					###2. consider different chrom with the 'site_chrm' and 'site_pos'???
					if algnmt.is_read1 == True:
						b_one_in_rep, b_both_in_rep = xalgnmt.is_TE_caused_discordant_pair(chrm, map_pos, mate_chrm, 
														mate_map_pos, xannotation)
					if b_one_in_rep:
						n_one_in_rep_region_pair += 1
					if b_both_in_rep:
						n_both_in_rep_region_pair += 1
					
					# m_new_candidate_sites[site_pos] = (
					#     str(n_barcode), str(n_discordant_pair), str(n_one_in_rep_region_pair), str(n_both_in_rep_region_pair),
					#     str(n_barcode_diff), str(n_barcode_share))
					m_new_candidate_sites[site_pos] = (str(n_barcode), str(n_discordant_pair), str(n_one_in_rep_region_pair), str(n_both_in_rep_region_pair))
		xbam.close_bam_file(barcode_bamfile)  ##close barcode bam file
		xbam.close_bam_file(bamfile)  ##close bam file
		
		##write out the combined results
		sf_candidate_list_disc = sf_candidate_list + global_values.DISC_SUFFIX_FILTER
		with open(sf_candidate_list_disc, "w") as fout_disc:
			for pos in m_new_candidate_sites:
				fout_disc.write(str(pos) + "\t")
				fout_disc.write(str(m_candidate_pos[pos]) + "\t")
				lth = len(m_new_candidate_sites[pos])
				for i in range(lth):
					fout_disc.write(str(m_new_candidate_sites[pos][i]) + "\t")
				fout_disc.write("\n")
				
	##filter out some false positive ones using the discordant reads
	def filter_candidate_sites_by_discordant_pairs(self, m_candidate_sites, iextend, i_is, f_dev, sf_annotation,
							n_one_in_rep_cutoff, n_both_in_rep_cutoff, sf_disc_working_folder):
		sf_disc_working_folder = self.working_folder + global_values.DISC_FOLDER
		if os.path.exists(sf_disc_working_folder) == False:
			cmd = "mkdir -p {0}".format(sf_disc_working_folder) # YW 2021/04/27 added -p
	    		#Popen(cmd, shell=True, stdout=PIPE).communicate()
			self.cmd_runner.run_cmd_small_output(cmd)
		sf_disc_working_folder += '/'
		self.output_candidate_sites_by_chrm(m_candidate_sites, sf_disc_working_folder, global_values.DISC_SUFFIX)
		l_chrm_records = []
		for chrm in m_candidate_sites:
			if len(chrm) > 5:  ###filter out those contigs!!!!!!! It's better to have a blacklist!!!!!!!!!!!!!!!!!!!!!!
				continue
			record = (chrm, self.sf_bam, self.sf_barcode_bam, iextend, i_is, f_dev, sf_annotation, sf_disc_working_folder, global_values.DISC_SUFFIX)
			l_chrm_records.append(record)
			# self.run_filter_by_discordant_pair_by_chrom(record) ##########################################################
		pool = Pool(self.n_jobs)
		pool.map(unwrap_self_filter_by_discordant, list(zip([self] * len(l_chrm_records), l_chrm_records)), 1)
		pool.close()
		pool.join()
		
		m_new_candidate_sites = {}
		for chrm in m_candidate_sites:  ####candidate site chromosome # read in by chrm
			sf_candidate_list_disc = sf_disc_working_folder + chrm + global_values.DISC_SUFFIX + global_values.DISC_SUFFIX_FILTER
			if os.path.exists(sf_candidate_list_disc) == False:
				continue
			with open(sf_candidate_list_disc) as fin_disc:
				for line in fin_disc:
					fields = line.split()
					pos = int(fields[0])
					n_barcode = int(fields[-6])
					n_discordant_pair = int(fields[-5])
					n_one_in_rep_region_pair = int(fields[-4])
					n_both_in_rep_region_pair = int(fields[-3])
					n_barcode_diff = int(fields[-2])
					n_barcode_share = int(fields[-1])
					
					if n_one_in_rep_region_pair < n_one_in_rep_cutoff:
						continue
					if n_both_in_rep_region_pair < n_both_in_rep_cutoff:
						continue
					
					if chrm not in m_new_candidate_sites:
						m_new_candidate_sites[chrm] = {}
					if pos not in m_new_candidate_sites[chrm]:
						if (chrm not in m_candidate_sites) or (pos not in m_candidate_sites[chrm]):
							continue
					n_clip = m_candidate_sites[chrm][pos][0]
					m_new_candidate_sites[chrm][pos] = (n_clip, n_barcode, n_discordant_pair, n_both_in_rep_region_pair, n_barcode_diff, n_barcode_share)
		return m_new_candidate_sites

	def run_filter_by_barcode_coverage(self, record):
		site_chrm1 = record[0]
		sf_bam = record[1]
		sf_barcode_bam = record[2]
		iextend = int(record[3])  ###extend some region on both sides in order to collect all barcodes
		i_cov_cutoff = int(record[4])
		sf_disc_working_folder = record[5]
		s_suffix = record[6]
		
		sf_candidate_list = sf_disc_working_folder + site_chrm1 + s_suffix
		if os.path.exists(sf_candidate_list) == False:
			return
		m_candidate_pos = {}
		with open(sf_candidate_list) as fin_list:
			for line in fin_list:
				fields = line.split()
				pos = int(fields[1])
				m_candidate_pos[pos] = "\t".join(fields[2:])
		
		bam_info = BamInfo(sf_bam, self.sf_reference)
		b_with_chr = bam_info.is_chrm_contain_chr()
		m_chrms = bam_info.get_all_reference_names()
		site_chrm = bam_info.process_chrm_name(site_chrm1, b_with_chr)
		if site_chrm not in m_chrms:
			return
		
		m_new_candidate_sites = {}
		xbam = XBamInfo(sf_bam, sf_barcode_bam, self.sf_reference)
		xbam.index_reference_name_id()
		bamfile = xbam.open_bam_file(sf_bam)  ##open bam file
		for site_pos in m_candidate_pos:  ####candidate site position # structure: {barcode:[alignmts]}
			if site_pos < iextend:
				continue
			set_barcodes = xbam.parse_barcodes_for_one_site(bamfile, site_chrm, site_pos, iextend)
			n_barcode = len(set_barcodes)
			if n_barcode > i_cov_cutoff: #if the barcode coverage is too high, then filter out the sites
				continue
			m_new_candidate_sites[site_pos] = n_barcode
		xbam.close_bam_file(bamfile)  ##close bam file
		
		##write out the combined results
		sf_candidate_list_disc = sf_candidate_list + global_values.DISC_SUFFIX_FILTER
		with open(sf_candidate_list_disc, "w") as fout_disc:
			for pos in m_new_candidate_sites:
				fout_disc.write(str(pos) + "\t")
				fout_disc.write(str(m_new_candidate_sites[pos]) + "\n")
	
	def filter_candidate_sites_by_barcode_coverage(self, m_candidate_sites, iextend, i_cov_cutoff):
		sf_disc_working_folder = self.working_folder + global_values.DISC_FOLDER
		if os.path.exists(sf_disc_working_folder) == False:
			cmd = "mkdir {0}".format(sf_disc_working_folder)
			self.cmd_runner.run_cmd_small_output(cmd)
		sf_disc_working_folder += '/'
		self.output_candidate_sites_by_chrm(m_candidate_sites, sf_disc_working_folder, global_values.DISC_SUFFIX)
		
		l_chrm_records = []
		for chrm in m_candidate_sites:
			if len(chrm) > 5:  ###filter out those contigs!!!!!!! It's better to have a blacklist!!!!!!!!!!!!!!!!!!!!!!
				continue
			record = (chrm, self.sf_bam, self.sf_barcode_bam, iextend, i_cov_cutoff, sf_disc_working_folder, global_values.DISC_SUFFIX)
			l_chrm_records.append(record)

		pool = Pool(self.n_jobs)
		pool.map(unwrap_self_filter_by_barcode_coverage, list(zip([self] * len(l_chrm_records), l_chrm_records)), 1)
		pool.close()
		pool.join()

		m_new_candidate_sites = {}
		for chrm in m_candidate_sites:  ####candidate site chromosome # read in by chrm
			sf_candidate_list_disc = sf_disc_working_folder + chrm + global_values.DISC_SUFFIX + global_values.DISC_SUFFIX_FILTER
			if os.path.exists(sf_candidate_list_disc) == False:
				continue
			with open(sf_candidate_list_disc) as fin_disc:
				for line in fin_disc:
					fields = line.split()
					pos = int(fields[0])
					n_barcode = int(fields[-1])
					
					if chrm not in m_new_candidate_sites:
						m_new_candidate_sites[chrm] = {}
					if pos not in m_new_candidate_sites[chrm]:
						if (chrm not in m_candidate_sites) or (pos not in m_candidate_sites[chrm]):
							continue
					n_clip = m_candidate_sites[chrm][pos][0]
					m_new_candidate_sites[chrm][pos] = (n_clip, n_barcode)
		return m_new_candidate_sites

	# YW 2020/08/03 github update: added everything related to sf_candidate_list_raw_disc
	# YW 2021/04/01 substitute sf_annotation w/ sf_annotation_Alu, sf_annotation_L1, sf_annotation_SVA
	# YW 2021/04/16 take out filtering by ratio of clustered disc reads
	def run_filter_by_discordant_pair_by_chrom_non_barcode(self, record):
		site_chrm1 = record[0]
		sf_bam = record[1]
		iextend = int(record[2])###extend some region on both sides in order to collect all barcodes
		i_is = int(record[3])
		f_dev = int(record[4])
		sf_annotation_Alu = record[5]
		sf_annotation_L1 = record[6]
		sf_annotation_SVA = record[7]
		sf_disc_working_folder = record[8]
		s_suffix = record[9]
		tmp_pos_folder = record[10]
		
		sf_candidate_list = sf_disc_working_folder + site_chrm1 + s_suffix
		sf_disc_pos = tmp_pos_folder + site_chrm1 + global_values.DISC_POS_SUFFIX # YW 2021/04/27 added to write chr level disc_pos
		if os.path.exists(sf_candidate_list) == False:
			return
		m_candidate_pos = {}
		with open(sf_candidate_list) as fin_list:
			for line in fin_list:
				fields = line.split()
				pos = int(fields[1])
				m_candidate_pos[pos] = "\t".join(fields[2:])
		os.remove(sf_candidate_list) # YW 2021/04/07 added to remove unused file
		
		bam_info = BamInfo(sf_bam, self.sf_reference)
		b_with_chr = bam_info.is_chrm_contain_chr()  # indicate whether the bam chrom has "chr" or not
		m_chrms = bam_info.get_all_reference_names()
		site_chrm = bam_info.process_chrm_name(site_chrm1, b_with_chr)
		if site_chrm not in m_chrms:
			return
		# print(site_chrm) ##########################################################################################
		xannotation_Alu = XAnnotation(sf_annotation_Alu)
		xannotation_Alu.set_with_chr(b_with_chr)
		xannotation_L1 = XAnnotation(sf_annotation_L1)
		xannotation_L1.set_with_chr(b_with_chr)
		xannotation_SVA = XAnnotation(sf_annotation_SVA)
		xannotation_SVA.set_with_chr(b_with_chr)
		# xannotation.load_rmsk_annotation()
		# xannotation.index_rmsk_annotation()
		i_min_copy_len=0
		boundary_extnd=0
		####For SVA, because the rmsk annotation is bad, so we extend a little bit
		# if global_values.IS_CALL_SVA is True:
		#     boundary_extnd=global_values.SVA_ANNOTATION_EXTND
		xannotation_Alu.load_rmsk_annotation_with_extnd_with_lenth_cutoff(boundary_extnd, i_min_copy_len)
		xannotation_Alu.index_rmsk_annotation_interval_tree()
		xannotation_L1.load_rmsk_annotation_with_extnd_with_lenth_cutoff(boundary_extnd, i_min_copy_len)
		xannotation_L1.index_rmsk_annotation_interval_tree()
		xannotation_SVA.load_rmsk_annotation_with_extnd_with_lenth_cutoff(global_values.SVA_ANNOTATION_EXTND, i_min_copy_len)
		xannotation_SVA.index_rmsk_annotation_interval_tree()
		
		bamfile = pysam.AlignmentFile(sf_bam, "rb", reference_filename=self.sf_reference)
		m_new_candidate_sites = {}
		m_raw_disc_sites={} # YW 2020/08/03 github update added
		m_chrm_ids=self._get_chrm_id_name(bamfile)
		for site_pos in m_candidate_pos:  ####candidate site position # structure: {barcode:[alignmts]}
			if site_pos < iextend:
				continue
			
			n_lraw_disc, n_left_discdt_Alu, n_left_discdt_L1, n_left_discdt_SVA, l_cluster = bam_info.cnt_collect_discordant_pairs(\
											bamfile, m_chrm_ids, site_chrm,
											site_pos - iextend, site_pos, i_is, f_dev, xannotation_Alu, xannotation_L1, xannotation_SVA, sf_disc_pos, "left")
			
			n_rraw_disc, n_right_discdt_Alu, n_right_discdt_L1, n_right_discdt_SVA, r_cluster = bam_info.cnt_collect_discordant_pairs(\
											bamfile, m_chrm_ids, site_chrm, site_pos + 1,
											site_pos + iextend, i_is, f_dev, xannotation_Alu, xannotation_L1, xannotation_SVA, sf_disc_pos, "right")
			m_new_candidate_sites[site_pos] = list(map(str, [n_left_discdt_Alu, n_right_discdt_Alu,
									 n_left_discdt_L1, n_right_discdt_L1,
									 n_left_discdt_SVA, n_right_discdt_SVA]))
			#if l_cluster[0]==True or r_cluster[0]==True:#for transductions, require at least one side form cluster
			# s_lcluster="0"
			r_lcluster = "0"
			s_lc_chrm = "-1"
			s_lc_pos = "-1"
			# if l_cluster[0]==True:
			#     s_lcluster="1"
			#     s_lc_chrm=l_cluster[1]
			#     s_lc_pos=str(l_cluster[2])
			r_lcluster = str(l_cluster[0])
			s_lc_chrm=l_cluster[1]
			s_lc_pos=str(l_cluster[2])
			# s_rcluster = "0"
			r_rcluster = "0"
			s_rc_chrm = "-1"
			s_rc_pos = "-1"
			# if r_cluster[0] == True:
			#     s_rcluster = "1"
			#     s_rc_chrm = r_cluster[1]
			#     s_rc_pos = str(r_cluster[2])
			r_rcluster = str(r_cluster[0])
			s_rc_chrm = r_cluster[1]
			s_rc_pos = str(r_cluster[2])
			
			m_raw_disc_sites[site_pos]=[str(n_lraw_disc), str(n_rraw_disc),
						    str(n_left_discdt_Alu), str(n_right_discdt_Alu),
						    str(n_left_discdt_L1), str(n_right_discdt_L1),
						    str(n_left_discdt_SVA), str(n_right_discdt_SVA),
						    r_lcluster, s_lc_chrm, s_lc_pos, r_rcluster, s_rc_chrm, s_rc_pos]
		bamfile.close()
		# YW 2021/04/26 added to remove empty files
		# if os.path.isfile(sf_disc_pos):
		# 	if os.path.getsize(sf_disc_pos)==0:
		# 		os.remove(sf_disc_pos)
####
		##write out the combined results
		sf_candidate_list_disc = sf_candidate_list + global_values.DISC_SUFFIX_FILTER
		with open(sf_candidate_list_disc, "w") as fout_disc:
			for pos in m_new_candidate_sites:
				fout_disc.write(str(pos) + "\t")
				# YW 2021/04/16 rewrote the following
				fout_disc.write("\t".join(m_new_candidate_sites[pos]) + "\n")
				# lth = len(m_new_candidate_sites[pos])
				# for i in range(lth):
				#     fout_disc.write(str(m_new_candidate_sites[pos][i]) + "\t")
				# fout_disc.write("\n")
		
		#save the raw disc results
		sf_candidate_list_raw_disc=sf_candidate_list + global_values.RAW_DISC_SUFFIX_FILTER
		with open(sf_candidate_list_raw_disc, "w") as fout_disc_raw:
			for pos in m_raw_disc_sites:
				fout_disc_raw.write(str(pos) + "\t")
				# YW 2021/04/16 rewrote the following
				fout_disc_raw.write("\t".join(m_raw_disc_sites[pos]) + "\n")
				# lth = len(m_raw_disc_sites[pos])
				# for i in range(lth):
				#     fout_disc_raw.write(str(m_raw_disc_sites[pos][i]) + "\t")
				# fout_disc_raw.write("\n")
	
	# YW 2020/08/03 github update: added everything related to sf_candidate_list_raw_disc
	# YW 2021/04/01 substitute sf_annotation w/ sf_annotation_Alu, sf_annotation_L1, sf_annotation_SVA
	# YW 2021/04/16 take out filtering by ratio of clustered disc reads
	# YW 2021/04/29 a shorter version of run_filter_by_discordant_pair_by_chrom_non_barcode: take out code to write sf_candidate_list + global_values.DISC_SUFFIX_FILTER
	def run_filter_by_discordant_pair_by_chrom_non_barcode_short(self, record):
		site_chrm1 = record[0]
		# print(f"run_filter_by_discordant_pair_by_chrom_non_barcode_short on {site_chrm1}!!!")
		sf_bam = record[1]
		iextend = int(record[2])###extend some region on both sides in order to collect all barcodes
		i_is = int(record[3])
		f_dev = int(record[4])
		sf_annotation_Alu = record[5]
		sf_annotation_L1 = record[6]
		sf_annotation_SVA = record[7]
		sf_disc_working_folder = record[8]
		s_suffix = record[9]
		tmp_pos_folder = record[10]
		
		sf_candidate_list = sf_disc_working_folder + site_chrm1 + s_suffix
		sf_disc_pos = tmp_pos_folder + site_chrm1 + global_values.DISC_POS_SUFFIX # YW 2021/04/27 added to write chr level disc_pos
		if os.path.exists(sf_candidate_list) == False:
			return
		m_candidate_pos = {}
		with open(sf_candidate_list) as fin_list:
			for line in fin_list:
				fields = line.split()
				pos = int(fields[1])
				m_candidate_pos[pos] = "\t".join(fields[2:])
		os.remove(sf_candidate_list) # YW 2021/04/07 added to remove unused file
		
		bam_info = BamInfo(sf_bam, self.sf_reference)
		b_with_chr = bam_info.is_chrm_contain_chr()  # indicate whether the bam chrom has "chr" or not
		m_chrms = bam_info.get_all_reference_names()
		site_chrm = bam_info.process_chrm_name(site_chrm1, b_with_chr)
		if site_chrm not in m_chrms:
			return
		# print(site_chrm) ##########################################################################################
		
		xannotation_Alu = XAnnotation(sf_annotation_Alu)
		xannotation_Alu.set_with_chr(b_with_chr)
		xannotation_L1 = XAnnotation(sf_annotation_L1)
		xannotation_L1.set_with_chr(b_with_chr)
		xannotation_SVA = XAnnotation(sf_annotation_SVA)
		xannotation_SVA.set_with_chr(b_with_chr)
		i_min_copy_len=0
		boundary_extnd=0
		xannotation_Alu.load_rmsk_annotation_with_extnd_with_lenth_cutoff(boundary_extnd, i_min_copy_len)
		xannotation_Alu.index_rmsk_annotation_interval_tree()
		xannotation_L1.load_rmsk_annotation_with_extnd_with_lenth_cutoff(boundary_extnd, i_min_copy_len)
		xannotation_L1.index_rmsk_annotation_interval_tree()
		xannotation_SVA.load_rmsk_annotation_with_extnd_with_lenth_cutoff(global_values.SVA_ANNOTATION_EXTND, i_min_copy_len)
		xannotation_SVA.index_rmsk_annotation_interval_tree()
		
		bamfile = pysam.AlignmentFile(sf_bam, "rb", reference_filename=self.sf_reference)
		m_new_candidate_sites = {}
		m_raw_disc_sites={} # YW 2020/08/03 github update added
		m_chrm_ids=self._get_chrm_id_name(bamfile)
		for site_pos in m_candidate_pos:  ####candidate site position # structure: {barcode:[alignmts]}
			if site_pos < iextend:
				continue

			n_lraw_disc, n_left_discdt_Alu, n_left_discdt_L1, n_left_discdt_SVA, l_cluster = bam_info.cnt_collect_discordant_pairs(\
										bamfile, m_chrm_ids, site_chrm,
										site_pos - iextend, site_pos, i_is, f_dev, xannotation_Alu, xannotation_L1, xannotation_SVA,sf_disc_pos, "left")
			
			n_rraw_disc, n_right_discdt_Alu, n_right_discdt_L1, n_right_discdt_SVA, r_cluster = bam_info.cnt_collect_discordant_pairs(\
										bamfile, m_chrm_ids, site_chrm, site_pos + 1,
										site_pos + iextend, i_is, f_dev, xannotation_Alu, xannotation_L1, xannotation_SVA, sf_disc_pos, "right")
			
			r_lcluster = "0"
			s_lc_chrm = "-1"
			s_lc_pos = "-1"
			r_lcluster = str(l_cluster[0])
			s_lc_chrm=l_cluster[1]
			s_lc_pos=str(l_cluster[2])
			
			r_rcluster = "0"
			s_rc_chrm = "-1"
			s_rc_pos = "-1"
			r_rcluster = str(r_cluster[0])
			s_rc_chrm = r_cluster[1]
			s_rc_pos = str(r_cluster[2])
			
			m_raw_disc_sites[site_pos]=[str(n_lraw_disc), str(n_rraw_disc),
						    str(n_left_discdt_Alu), str(n_right_discdt_Alu),
						    str(n_left_discdt_L1), str(n_right_discdt_L1),
						    str(n_left_discdt_SVA), str(n_right_discdt_SVA),
						    r_lcluster, s_lc_chrm, s_lc_pos, r_rcluster, s_rc_chrm, s_rc_pos]
		bamfile.close()
####
		##write out the combined results
		#save the raw disc results
		sf_candidate_list_raw_disc=sf_candidate_list + global_values.RAW_DISC_SUFFIX_FILTER
		with open(sf_candidate_list_raw_disc, "w") as fout_disc_raw:
			for pos in m_raw_disc_sites:
				fout_disc_raw.write(str(pos) + "\t")
				# YW 2021/04/16 rewrote the following
				fout_disc_raw.write("\t".join(m_raw_disc_sites[pos]) + "\n")
				

	###This one feed in the normal illumina data, and count the discordant pairs of the left and right regions
	# YW 2021/04/01 changed the function input
	# YW 2021/04/17 update dict elements: ndisc for Alu/L1/SVA
	# YW 2021/04/27 add code to aggregate disc_pos from different chrs
	def filter_candidate_sites_by_discordant_pairs_non_barcode(self, m_candidate_sites, iextend, i_is, f_dev,
															   sf_annotation_Alu, sf_annotation_L1, sf_annotation_SVA, n_discordant_cutoff, sf_disc_pos_tmp, tmp_pos_folder):
		#####################################
		# YW 2021/04/27 moved this to filter_candidate_sites_by_discordant_pairs_multi_alignmts
		# sf_disc_working_folder = self.working_folder + global_values.DISC_FOLDER
		# if os.path.exists(sf_disc_working_folder) == False:
		# 	cmd = "mkdir {0}".format(sf_disc_working_folder)
		# 	#Popen(cmd, shell=True, stdout=PIPE).communicate()
		# 	self.cmd_runner.run_cmd_small_output(cmd)
		# sf_disc_working_folder += '/'
		#####################################
		self.output_candidate_sites_by_chrm(m_candidate_sites, sf_disc_working_folder, global_values.DISC_SUFFIX)
		xchrm=XChromosome()
		l_chrm_records = []
		for chrm in m_candidate_sites:
			###filter out those contigs!!!!!!!
			if xchrm.is_decoy_contig_chrms(chrm)==True:
				continue
			
			l_chrm_records.append(
					(chrm, self.sf_bam, iextend, i_is, f_dev, sf_annotation_Alu, sf_annotation_L1, sf_annotation_SVA, sf_disc_working_folder, global_values.DISC_SUFFIX, tmp_pos_folder)) # YW 2021/04/27 added sf_disc_fa_tmp, tmp_pos_folder to save disc_pos
		
		pool = Pool(self.n_jobs)
		pool.map(unwrap_self_filter_by_discordant_non_barcode, list(zip([self] * len(l_chrm_records), l_chrm_records)), 1)
		pool.close()
		pool.join()
		
		##########################################################
		# YW 2021/04/27 added to aggregate disc pos from each chr (within the same sample)
		with open(sf_disc_pos_tmp, 'w') as fout_merged_disc_pos:
			for chrm in m_candidate_sites:
				if xchrm.is_decoy_contig_chrms(chrm)==True:
					continue
				sf_disc_pos = tmp_pos_folder + chrm + global_values.DISC_POS_SUFFIX
				if os.path.exists(sf_disc_pos):
					with open(sf_disc_pos) as fin_disc_pos:
						for line in fin_disc_pos:
							fout_merged_disc_pos.write(line)
					os.remove(sf_disc_pos)
		##########################################################
		
		m_new_candidate_sites = {}
		for chrm in m_candidate_sites:  ####candidate site chromosome # read in by chrm
			sf_candidate_list_disc = sf_disc_working_folder + chrm + global_values.DISC_SUFFIX + global_values.DISC_SUFFIX_FILTER
			if os.path.exists(sf_candidate_list_disc) == False:
				continue
			with open(sf_candidate_list_disc) as fin_disc:
				for line in fin_disc:
					fields = line.split()
					pos = int(fields[0])
					n_disc_left_Alu = int(fields[1])
					n_disc_right_Alu = int(fields[2])
					n_disc_left_L1 = int(fields[3])
					n_disc_right_L1 = int(fields[4])
					n_disc_left_SVA = int(fields[5])
					n_disc_right_SVA = int(fields[6])
					
					#	##Here require both left and right discordant pairs
					# if n_disc_left < n_discordant_cutoff and n_disc_right < n_discordant_cutoff:
					#     continue
					if (n_disc_left_Alu + n_disc_right_Alu < n_discordant_cutoff) and (n_disc_left_L1 + n_disc_right_L1 < n_discordant_cutoff) and (n_disc_left_SVA + n_disc_right_SVA < n_discordant_cutoff):
						continue
					
					if chrm not in m_new_candidate_sites:
						m_new_candidate_sites[chrm] = {}
					if pos not in m_new_candidate_sites[chrm]:
						# n_clip = m_candidate_sites[chrm][pos][0]
						# m_new_candidate_sites[chrm][pos] = (n_clip, n_disc_left, n_disc_right)
						m_new_candidate_sites[chrm][pos] = [n_disc_left_Alu, n_disc_right_Alu, n_disc_left_L1, n_disc_right_L1, n_disc_left_SVA, n_disc_right_SVA]
			os.remove(sf_candidate_list_disc) # YW 2021/04/07 added to remove merged file
		
		#RAW_DISC_SUFFIX_FILTER
		m_raw_disc_sites = {}  # for each site, save the raw discordant reads
		for chrm in m_candidate_sites:  ####candidate site chromosome # read in by chrm
			sf_candidate_list_raw_disc = sf_disc_working_folder + chrm + global_values.DISC_SUFFIX + global_values.RAW_DISC_SUFFIX_FILTER
			if os.path.exists(sf_candidate_list_raw_disc) == False:
				continue
			with open(sf_candidate_list_raw_disc) as fin_disc:
				for line in fin_disc:
					fields = line.split()
					pos = int(fields[0])
					n_raw_disc_left = int(fields[1])
					n_raw_disc_right = int(fields[2])
					n_disc_left_Alu = int(fields[3])
					n_disc_right_Alu = int(fields[4])
					n_disc_left_L1 = int(fields[5])
					n_disc_right_L1 = int(fields[6])
					n_disc_left_SVA = int(fields[7])
					n_disc_right_SVA = int(fields[8])
					r_lcluster=float(fields[9])
					s_lc_chrm=fields[10]
					s_lc_pos=fields[11]
					r_rcluster=float(fields[12])
					s_rc_chrm=fields[13]
					s_rc_pos=fields[14]
					
					if chrm not in m_raw_disc_sites:
						m_raw_disc_sites[chrm] = {}
						if pos not in m_raw_disc_sites[chrm]:
						# n_clip = m_candidate_sites[chrm][pos][0]
						# m_new_candidate_sites[chrm][pos] = (n_clip, n_disc_left, n_disc_right)
							m_raw_disc_sites[chrm][pos] = [n_raw_disc_left, n_raw_disc_right, n_disc_left_Alu, n_disc_right_Alu,
										       n_disc_left_L1, n_disc_right_L1, n_disc_left_SVA, n_disc_right_SVA,
										       r_lcluster, s_lc_chrm, s_lc_pos, r_rcluster, s_rc_chrm, s_rc_pos]
			os.remove(sf_candidate_list_raw_disc) # YW 2021/04/07 added to remove merged file
		return m_new_candidate_sites, m_raw_disc_sites
	
	# YW 2021/04/29 shorten the filter_candidate_sites_by_discordant_pairs_non_barcode function to get rid of m_new_candidate_sites output
	def filter_candidate_sites_by_discordant_pairs_non_barcode_short(self, m_candidate_sites, iextend, i_is, f_dev,
																	 sf_annotation_Alu, sf_annotation_L1, sf_annotation_SVA, n_discordant_cutoff, sf_disc_pos_tmp, sf_disc_working_folder, tmp_pos_folder):
		self.output_candidate_sites_by_chrm(m_candidate_sites, sf_disc_working_folder, global_values.DISC_SUFFIX)
		xchrm=XChromosome()
		l_chrm_records = []
		for chrm in m_candidate_sites:
			###filter out those contigs!!!!!!!
			if xchrm.is_decoy_contig_chrms(chrm)==True:
				continue
			
			l_chrm_records.append(
					(chrm, self.sf_bam, iextend, i_is, f_dev, sf_annotation_Alu, sf_annotation_L1, sf_annotation_SVA, sf_disc_working_folder, global_values.DISC_SUFFIX, tmp_pos_folder)) # YW 2021/04/27 added tmp_pos_folder to save disc_pos
		
		pool = Pool(self.n_jobs)
		pool.map(unwrap_self_filter_by_discordant_non_barcode_short, list(zip([self] * len(l_chrm_records), l_chrm_records)), 1)
		pool.close()
		pool.join()
		
		##########################################################
		# YW 2021/04/27 added to aggregate disc pos from each chr (within the same sample)
		with open(sf_disc_pos_tmp, 'w') as fout_merged_disc_pos:
			for chrm in m_candidate_sites:
				if xchrm.is_decoy_contig_chrms(chrm)==True:
					continue
				sf_disc_pos = tmp_pos_folder + chrm + global_values.DISC_POS_SUFFIX
				if os.path.exists(sf_disc_pos):
					with open(sf_disc_pos) as fin_disc_pos:
						for line in fin_disc_pos:
							fout_merged_disc_pos.write(line)
					os.remove(sf_disc_pos)
		##########################################################
		
		#RAW_DISC_SUFFIX_FILTER
		m_raw_disc_sites = {}  # for each site, save the raw discordant reads
		for chrm in m_candidate_sites:  ####candidate site chromosome # read in by chrm
			sf_candidate_list_raw_disc = sf_disc_working_folder + chrm + global_values.DISC_SUFFIX + global_values.RAW_DISC_SUFFIX_FILTER
			if os.path.exists(sf_candidate_list_raw_disc) == False:
				continue
			with open(sf_candidate_list_raw_disc) as fin_disc:
				for line in fin_disc:
					fields = line.split()
					pos = int(fields[0])
					n_raw_disc_left = int(fields[1])
					n_raw_disc_right = int(fields[2])
					n_disc_left_Alu = int(fields[3])
					n_disc_right_Alu = int(fields[4])
					n_disc_left_L1 = int(fields[5])
					n_disc_right_L1 = int(fields[6])
					n_disc_left_SVA = int(fields[7])
					n_disc_right_SVA = int(fields[8])
					r_lcluster=float(fields[9])
					s_lc_chrm=fields[10]
					s_lc_pos=fields[11]
					r_rcluster=float(fields[12])
					s_rc_chrm=fields[13]
					s_rc_pos=fields[14]
					
					if chrm not in m_raw_disc_sites:
						m_raw_disc_sites[chrm] = {}
					if pos not in m_raw_disc_sites[chrm]:
						m_raw_disc_sites[chrm][pos] = [n_raw_disc_left, n_raw_disc_right, n_disc_left_Alu, n_disc_right_Alu,
									       n_disc_left_L1, n_disc_right_L1, n_disc_left_SVA, n_disc_right_SVA,
									       r_lcluster, s_lc_chrm, s_lc_pos, r_rcluster, s_rc_chrm, s_rc_pos]
			os.remove(sf_candidate_list_raw_disc) # YW 2021/04/07 added to remove merged file
		return m_raw_disc_sites

	def output_candidate_sites_by_chrm(self, m_candidate_list, sf_folder, s_suffix):
		for chrm in m_candidate_list:
			with open(sf_folder + chrm + s_suffix, "w") as fout_chrm:
				for pos in m_candidate_list[chrm]:
					# lth = len(m_candidate_list[chrm][pos])
					# fout_chrm.write(chrm + "\t" + str(pos) + "\t")
					fout_chrm.write("\t".join([chrm, str(pos), ""]))
					fout_chrm.write("\t".join([str(i) for i in m_candidate_list[chrm][pos]]) + "\n")
					# for i in range(lth):
					#     fout_chrm.write(str(m_candidate_list[chrm][pos][i]) + "\t")
					# fout_chrm.write("\n")

	###output the candidate list in a file
	def output_candidate_sites(self, m_candidate_list, sf_out):
		with open(sf_out, "w") as fout_candidate_sites:
			for chrm in m_candidate_list:
				for pos in m_candidate_list[chrm]:
					lth = len(m_candidate_list[chrm][pos])
					fout_candidate_sites.write(chrm + "\t" + str(pos) + "\t")
					for i in range(lth):
						s_feature = str(m_candidate_list[chrm][pos][i])
						fout_candidate_sites.write(s_feature + "\t")
					fout_candidate_sites.write("\n")

	####Two situations will be kept:
	# 1. the site have both left and right clipped reads, and mate reads fall in repeat region
	# 2. The site has nearby sites, that form a cluster and the cluster have both left and right clpped reads
	# def is_candidate_clip_position(self, pos, m_clip_pos, iextend, n_left_cutoff, n_right_cutoff, n_mate_cutoff):
	#     if

	####according to the clipped reads (left, right) clip:
	# 1. the clip position either has left_clipped and right_clipped reads, or
	# 2. the clip position has left (right) clipped reads, and nearby clipped position has right (left) clip reads
	####According to barcode:
	# 1. at least 20 barcodes are shared between the left and right regions
	# 2.
	def first_stage_filter(self, sf_working_folder):
		bam_info = BamInfo(self.sf_bam, self.sf_reference)
		m_chrm_names = bam_info.get_all_reference_names()
		m_all_sites = {}
		
		for chrm in m_chrm_names:
			m_clip_pos = {}
			sf_clip_sites = sf_working_folder + "{0}{1}".format(chrm, global_values.CLIP_POS_SUFFIX)
			if os.path.exists(sf_clip_sites) == False:
				continue
			sf_disc_sites = sf_working_folder + global_values.DISC_FOLDER + "/" + "{0}{1}{2}".format(chrm, global_values.DISC_SUFFIX,
					global_values.DISC_SUFFIX_FILTER)
			if os.path.exists(sf_disc_sites) == False:
				continue
		with open(sf_clip_sites) as fin_clip:  # read in the sites from file
			for line in fin_clip:
				fields = line.split()
				pos = int(fields[0])
				n_left_clip = int(fields[1])
				n_right_clip = int(fields[2])
				n_mate_in_rep = int(fields[3])
				m_clip_pos[pos] = (n_left_clip, n_right_clip, n_mate_in_rep)
		
		##TDList: Need to process the clipped reads
		# need to collect the clipped part, and re-align to the reference
		
		# 1. find all the clipped position, and then re-align the clipped part ###count the number of clipped reads
		# 1.1 number of left clipped, and number of right clipped reads
		# 1.2 number of (discordant, clipped) reads
		# 2. for each clip position, find all the discordant pairs             ###count the number of discordant pairs
		# 3.
		
		# 4. Align the assembled contigs to the reference (align the repeat copies to the assembled contigs????)
		# 5. call out the events
		# 6.
		
		# 7.
		####

	# def merge_sites_features(self, s_working_folder, m_sites_barcode, sf_out):
	def _get_chrm_id_name(self, samfile):
		m_chrm = {}
		references = samfile.references
		for schrm in references:
			chrm_id = samfile.get_tid(schrm)
			m_chrm[chrm_id] = schrm
		m_chrm[-1] = "*"
		return m_chrm

####
    #import os, tempfile
    # def symlink(self, target, link_name, overwrite=False):
    #     '''
    #     Create a symbolic link named link_name pointing to target.
    #     If link_name exists then FileExistsError is raised, unless overwrite=True.
    #     When trying to overwrite a directory, IsADirectoryError is raised.
    #     '''
    #
    #     if not overwrite:
    #         os.symlink(target, linkname)
    #         return
    #
    #     # os.replace() may fail if files are on different filesystems
    #     link_dir = os.path.dirname(link_name)
    #
    #     # Create link to target with temporary filename
    #     while True:
    #         temp_link_name = tempfile.mktemp(dir=link_dir)
    #
    #         # os.* functions mimic as closely as possible system functions
    #         # The POSIX symlink() returns EEXIST if link_name already exists
    #         # https://pubs.opengroup.org/onlinepubs/9699919799/functions/symlink.html
    #         try:
    #             os.symlink(target, temp_link_name)
    #             break
    #         except FileExistsError:
    #             pass
    #
    #     # Replace link_name with temp_link_name
    #     try:
    #         # Pre-empt os.replace on a directory with a nicer message
    #         if os.path.isdir(link_name):
    #             raise IsADirectoryError(f"Cannot symlink over existing directory: '{link_name}'")
    #         os.replace(temp_link_name, link_name)
    #     except:
    #         if os.path.islink(temp_link_name):
    #             os.remove(temp_link_name)
    #         raise
