#!/bin/python3
# 2021/05/26
##@@author: Yilan (Elain) Wang, Harvard University
##@@contact: yilanwang@g.harvard.edu

# REMEMBER TO CHANGE THE EMAIL ADDRESS IN the run_sbatch_scripts function IF YOU USE THIS SCRIPT!!!!!

'''
To speed up the xTea-ML feature extraction process for training, validation, and testing
Generate python, bash, and sbatch scripts to run realignment of clipped/discordant reads to repeat cns in parallel
This module is run for each bam file (runs multiple times if there are more than one bam files for a sample)
'''
# import python standard libraries
import os
import sys

# import modules (source_scripts is the symbolic link to the xTea-ML/xtea directory)
import global_values
from bwa_align import *
from cmd_runner import *

class CNS_PARALLEL():
	
	def __init__(self, n_jobs, working_folder, sf_fq, idx_bam):
		self.n_jobs = n_jobs
		self.working_folder = working_folder
		self.sf_fq = sf_fq
		self.idx_bam = idx_bam
		self.cmd_runner = CMD_RUNNER()
	
	def gnrt_clip_py_scripts(self, rep_type, sf_rep_cns, sf_rep, sf_algnmt):
		py_scr = self.working_folder + rep_type + "_clip_cns_realn.py"
		done_f = self.working_folder + rep_type + "_clip_cns_realn.done"
		fail_f = self.working_folder + rep_type + "_clip_cns_realn.fail"
		if os.path.exists(done_f):
			os.remove(done_f)
		if os.path.exists(fail_f):
			os.remove(fail_f)
		with open(py_scr, 'w') as fout:
			scmd = "#!/bin/python3\n\n"
			scmd += "import sys\n"
			scmd += f"sys.path.insert(1, \"{self.working_folder}/source_scripts\")\n" # to resolve the issues with importing within bwa_align.py
			scmd += "import global_values\n"
			scmd += "from bwa_align import *\n\n"
			# scmd += "import source_scripts.global_values\n"
			# scmd += "from source_scripts.bwa_align import *\n"
			scmd += f"bwa_align = BWAlign(global_values.BWA_PATH, global_values.BWA_REALIGN_CUTOFF, {self.n_jobs})\n"
			scmd += f"bwa_align.two_stage_realign(\"{sf_rep_cns}\", \"{sf_rep}\", \"{self.sf_fq}\", \"{sf_algnmt}\")\n"
			fout.write(scmd)
		return py_scr, done_f, fail_f
	
	def gnrt_disc_py_scripts(self, rep_type, sf_rep_cns, sf_algnmt):
		py_scr = self.working_folder + rep_type + "_disc_cns_realn.py"
		done_f = self.working_folder + rep_type + "_disc_cns_realn.done"
		fail_f = self.working_folder + rep_type + "_disc_cns_realn.fail"
		if os.path.exists(done_f):
			os.remove(done_f)
		if os.path.exists(fail_f):
			os.remove(fail_f)
		with open(py_scr, 'w') as fout:
			scmd = "#!/bin/python3\n\n"
			scmd += "import sys\n"
			scmd += f"sys.path.insert(1, \"{self.working_folder}/source_scripts\")\n" # to resolve the issues with importing within bwa_align.py
			scmd += "import global_values\n"
			scmd += "from bwa_align import *\n\n"
			# scmd += "import source_scripts.global_values\n"
			# scmd += "from source_scripts.bwa_align import *\n"
			scmd += f"bwa_align = BWAlign(global_values.BWA_PATH, global_values.BWA_REALIGN_CUTOFF, {self.n_jobs})\n"
			scmd += f"bwa_align.realign_disc_reads(\"{sf_rep_cns}\", \"{self.sf_fq}\", \"{sf_algnmt}\")\n"
			fout.write(scmd)
		return py_scr, done_f, fail_f
	
	def run_sbatch_scripts(self, scr_list, option):
		if option != "clip" and option != "disc":
			raise NotImplementedError
		bash_script = self.working_folder + option + "_cns_realn.sh"
		num = len(scr_list)
		with open(bash_script, 'w') as fout:
			fout.write("#!/bin/bash\n\n")
			for scr in scr_list:
				py_scr, done_f, fail_f = scr
				rep_type = os.path.basename(py_scr).split("_")[0]
				sbatch_script = self.working_folder + rep_type + "_" + option + "_cns_realn.sh"
				with open(sbatch_script, 'w') as sub_fout:
					sheader = "#!/bin/bash\n\n"
					if option == "clip":
						sheader += f"#SBATCH -p {global_values.C_REALIGN_PARTITION}\n"
						sheader += f"#SBATCH -t {global_values.C_REALIGN_TIME}\n"
						sheader += f"#SBATCH --mem={global_values.C_REALIGN_MEMORY}G\n"
					else:
						sheader += f"#SBATCH -p {global_values.D_REALIGN_PARTITION}\n"
						sheader += f"#SBATCH -t {global_values.D_REALIGN_TIME}\n"
						sheader += f"#SBATCH --mem={global_values.D_REALIGN_MEMORY}G\n"
					sheader += "#SBATCH -N 1\n"
					sheader += f"#SBATCH -c {self.n_jobs}\n"
					sheader += f"#SBATCH -J {self.idx_bam}_{rep_type}_{option}_cns_realn\n"
					sheader += f"#SBATCH -o {self.working_folder}/{self.idx_bam}_{rep_type}_{option}_cns_%j.out\n"
					sheader += "#SBATCH --mail-type=FAIL\n"
					sheader += "#SBATCH --mail-user=yilanwang@g.harvard.edu\n\n"
					scmd = f"python3 {py_scr}\n"
					scmd += f"STATUS=$?\n"
					scmd += f"if [ \"$STATUS\" -eq 0 ]; then touch {done_f}; else touch {fail_f}; fi\n"
					sub_fout.write(sheader)
					sub_fout.write(scmd)
				fout.write(f"sbatch {sbatch_script}\n")
		bash_cmd = f"bash {bash_script}"
		self.cmd_runner.run_cmd_small_output(bash_cmd)
		print(f"Successfully launched {num} {option} cns realignment job(s)!")
