#!/bin/python3
# 2021/07/27
##@@author: Yilan (Elain) Wang, Harvard University
##@@contact: yilanwang@g.harvard.edu
'''
1. modified from 3.0coor_lift.py, put everything into a python class, used dynamic programming(?) to accelerate the runs
2. convert TE coordinates from ctrl bam to test bam
Note: this is not an exact inverse function of 3.0coor_lift.py parse_result_ML because insertions located within gold std set TEs are processed differently
3. subtract TE coordinates in test bam overlapping those in ctrl bam after coordinate lifting
'''
# import python standard libraries
import os
import sys

# import modules
from cmd_runner import *

class Coor_Lift():
	
	def __init__(self, input, output, ref, error_margin=15):
		self.input = input
		self.output = output
		self.ref = ref
		self.error_margin = error_margin
		self.cmd_runner = CMD_RUNNER()
	
	
	@staticmethod
	def write_output(f_out, *argv):
		f_out.write("\t".join(argv)+'\n')
	
	
	@staticmethod
	def load_ref(ref):
		'''
		load ref bed into a dict for faster processing
		now compatible when there are different chr present in ref_bed
		'''
		ref_dict = {}
		with open(ref, 'r') as gold_std:
			for LINES in gold_std:
				LINE = LINES.rstrip().split('\t')
				CHR = LINE[0]
				START = int(LINE[1])
				END = int(LINE[2])
				INTERVAL = int(LINE[4])
				if CHR not in ref_dict:
					ref_dict[CHR] = [{'START':START, 'END':END, 'INTERVAL':INTERVAL}]
				else:
					ref_dict[CHR].append({'START':START, 'END':END, 'INTERVAL':INTERVAL})
		return ref_dict
			
	
	def run_coor_lift(self, mode="ctrl"):
		'''
		convert TE coordinates from ctrl bam to test bam
		'''
		if mode != "ctrl":
			raise NotImplementedError
		
		ref_dict = self.load_ref(self.ref)
		chr_list = list(ref_dict.keys())
		ref_chr_pos = dict.fromkeys(chr_list, 0) # remember which line number we have reached for each chr (avoid starting from scratch)
		to_change_sum = dict.fromkeys(chr_list, 0) # remember how much to add/subtract based on sum of INTERVAL (avoid starting from scratch)
		print(f"Performing coordinate lifting on {mode} and chrm:")
		print(chr_list)
		# make sure the input file is sorted
		self.cmd_runner.run_cmd_to_file(f"sort -V -k 1,2 {self.input}", self.input + ".sorted")
		with open(self.input + ".sorted", 'r') as to_lift, open(self.output, 'w') as f_out:
			for lines in to_lift:
				if lines[0] == "#":
					continue
				line = lines.rstrip().split('\t')
				chr = line[0]
				s = int(line[1])
				if chr not in ref_dict: # may need to consider processing chr names
					self.write_output(f_out, chr, str(s), str(s+1))
				else: # for chr that contains gold std set
					if ref_chr_pos[chr] < len(ref_dict[chr]): # haven't iterated through all the gold std set pos in chr
						cnt = 0 # count how many more gold std set we have gone through
						for gold_std_set in ref_dict[chr][ref_chr_pos[chr]:]: # start iterating from the last position in the previous
							START = gold_std_set['START']
							END = gold_std_set['END']
							INTERVAL = gold_std_set['INTERVAL']
							if mode == "ctrl":
								if s > START - self.error_margin: # for now, do not subtract the current INTERVAL for those landing w/in gold std set
									if s <= END + self.error_margin: # may need to check again
										self.write_output(f_out, chr, str(s-to_change_sum[chr]), str(s+1-to_change_sum[chr]))
										break
									to_change_sum[chr] += INTERVAL
									cnt += 1
								else:
									self.write_output(f_out, chr, str(s-to_change_sum[chr]), str(s+1-to_change_sum[chr]))
									break
						ref_chr_pos[chr] += cnt # update the position of gold std set we should start with
						
					else: # have iterated through all the gold std set pos in chr
						if mode == 'ctrl':
							self.write_output(f_out, chr, str(s-to_change_sum[chr]), str(s+1-to_change_sum[chr]))
		# clean up the intermediate file
		os.remove(self.input + ".sorted")
	
	
	def sort_subtract_overlap(self, ctrl_file):
		'''
		subtract TE coordinates in test bam overlapping those in ctrl bam after coordinate lifting
		test bam only
		'''
		print(f"Removing overlapping coordinate in {self.input} from {ctrl_file}...")
		if not os.path.isfile(ctrl_file):
			sys.exit(f"{ctrl_file} does not exist. Please run xTea-ML on the corresponding control bam file first!")
		sorted_input = self.input + "_sorted" # do not use .sorted as this is self.output
		sorted_ctrl = ctrl_file + ".sorted"
		# make sure the input file is sorted
		self.cmd_runner.run_cmd_to_file(f"sort -V -k 1,2 {self.input}", sorted_input)
		self.cmd_runner.run_cmd_to_file(f"sort -V -k 1,2 {ctrl_file}", sorted_ctrl)
		self.cmd_runner.run_cmd_to_file("awk -F'\\t' '{OFS=\"\\t\"; print $1,$2,$2+1,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24}' " + f"{sorted_input} | intersectBed -v -sorted -wa -a stdin -b {sorted_ctrl}", self.output)
		# clean up intermediate files
		os.remove(sorted_input)
		os.remove(sorted_ctrl)
						