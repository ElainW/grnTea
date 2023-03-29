##11/27/2017
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu

'''
Upgrade 11/04/2018
1. Put all the pubic shared variables to a single file global_variables.py
2. For clip and disc have two different mapping-quality cutoff now
3. Keep the short poly-A clipped reads in the "clip" step
4. Change the align clipped reads to repeat copies to "two stages" alignment
5. Popen.communicate() has a potential bug when the "stdout" file is large. Now change it redirect to a file,
    instead of stdout
6. Add the gene annotation information to each candidate
####11/12/2018 --solved on 11/14/2018
Bug need to fix: when have several bams for one individual, the pub_folder link may have the same file name,
this is not correct.
####

####@@CC 01/03/2018
Todo List: 1. clip step requires lots of memory. 150G for 10 cores. This should be improved.
Update: 01/09/2018: The issue is caused at the "ClipReadInfo::cnt_clip_part_aligned_to_rep()" function.
Solution:
        1)Skip the load all into memory step, instead only load those has "left(right) realigned" ones. This is
            the current version used.
        2)For the re-aligned sam file, parse out the chrm, positions, and write to a file, and sort the file,
          then, count the number of supported reads. This version should be more memory efficient, but hasn't been
          fully tested.

####
09/10/2018
Todo List:
1. For orphan (or partial) transduction events:
    1) In "clip" step, change "mate-in-rep" to "mate-is-discordant"
    3) In "disc" step, change "mate-in-rep" to "mate-is-discordant", then split to two outputs:
3. output to vcf file


10/10/2018
Trio checking module:
Need to check whether the clipped parts are same or not!

11/06/2018:
Add a coverage estimation module. (need to sum up the coverage of different bams) --To do list
Set an automatic cutoff setting module that associated with read depth, e.g. (solved 01/05/2019)
30-40X: clip:3, disc:5
50-60: clip:4, disc:6
70-100: clip:5, disc:8

Also, need to set the maximum coverage cutoff (current one in filtering step, default set as 8*30=240X)
Similarly, check_disc_sample_consistency(), should use an arrary of cutoff value.

11/20/2018 Two major modules now affact the speed a lot, especiallly when there lots of clipped candidate sites (for 10X):
1. Align the clipped reads to repeat copies. (Already have been improved, but still lots of rooms)
2. Check discordant reads (second step) step, current module is split by chrm, and still have room to improve!!!

To do list:
12/11/2018 (solved)
Check the phred score for the clipped parts, and if the score is low, then discard the clipped part
#Q: how to set the cutoff for this value? By default set as 15

12/12/2018 (solved on 12/31/18)
Collect all the background information
#nearby un-related clipped reads
#nearby un-related or low mapping quality discordant reads
#fully mapped reads at the breakpoint
#concordant reads near the breakpoint
##reads contain small indels or lots of mismatches.... (how to calc this) ???
##
12/26/2018 (solved on 12/31/18)
Revise: polyA part now changed to only check sub-clipped reads

##
12/31/18
Add one option for cleaning tmp files.

01/31/2019
Potential issue:
1. When reads are long say: 250bp, but the insertion is short like 200bp, then there will be no/small-number of discordant reads.
2. This need to be solved when adding the long reads module.

02/28/2019
Improvement:
1. Define high confident regions (discordant reads aligned to these regions will not be counted)
2. Add repeatmasker annotation in the final output
3. If one site is called as both Alu and SVA (high confident), then select as SVA
4. When checking whether neighor region is coverage island or not, count all reads, not only those with high mapq.
5. Check the TSD assembly? To do further filtering.
    For those have two side clipped reads support, left/right breakpoint is supported by left/right clipped reads only
6. Check the standard deviation for breakpoint clipped reads (maybe also the discordant reads?), especially for 10X data

03/24/2019
Improvement:
For cram format, if we don't need to parse quality, seq fields, then some species process can improve the speed:
See here: https://brentp.github.io/post/cram-speed/
bam_list.append(pysam.AlignmentFile(b,
                mode='rc',reference_filename=ref_fasta,format_options=["required_fields=7167"]))
https://gist.github.com/brentp/213f214dd5677dd447150e62e2e360f5

04/17/19 (done, need to collect further information )
Collect the left-disc-rc/non-rc, and right-dis-rc/non-rc features and save them in cns file;
Note: these are from aligned ones, also need to check those unaligned ones
####
Also need to collect the: direction of the background ones !!!!!
####

####
04/17/19 Add the post filtering module (done)
####
'''
####

####
import os
from shutil import copyfile
import argparse
import subprocess # to do line count of the disc read output
from glob import glob # to get a list of all the tmp disc read output files

import global_values
from x_TEI_locator import *
from x_local_assembly import *
from x_intermediate_sites import *
from x_reference import *
from x_clip_disc_filter import *
from x_genotype_feature import *
from x_basic_info import *
from x_parameter import *
# from x_post_filter import *
# from x_joint_calling import *
# from x_igv import *
# from x_genotype_classify import * # new
from extract_features import * # YW 2021/05/10 added this
from coor_lift import * # YW 2021/07/27 added this
from cmd_runner import * # YW 2021/10/27 added this


def parse_arguments():
    parser = argparse.ArgumentParser("Main script of xTea_ML")
    parser.add_argument("-C", "--clip",
                        action="store_true", dest="clip", default=False,
                        help="Call candidate TEI sites from clipped reads")
    parser.add_argument("-S", "--single",
                        action="store_true", dest="single", default=False,
                        help="Call clip positions from single-end reads")
    parser.add_argument("-D", "--discordant",
                        action="store_true", dest="discordant", default=False,
                        help="Filter with discordant paired end reads")
    parser.add_argument("-N", "--filter_csn",
                        action="store_true", dest="filter_csn", default=False,
                        help="Filter out candidate sites from map position on consensus")
    parser.add_argument("--resume",
                        action="store_true", dest="resume", default=False,
                        help="Resume the running, which will skip the step if output file already exists!")
    parser.add_argument("--mit",
                        action="store_true", dest="mit", default=False,
                        help="Indicate call mitochondrion insertion")
    parser.add_argument("--dna",
                        action="store_true", dest="dna", default=False,
                        help="Not RNA mediated insertion (no polyA)")
    parser.add_argument("--cbs",
                        action="store_true", dest="cbs", default=False,
                        help="check by sample")#whether check by sample
    parser.add_argument("--sva",
                        action="store_true", dest="sva", default=False,
                        help="For SVA insertion calling")
    # YW 2020/08/16 added this option
    parser.add_argument("--l1",
                        action="store_true", dest="l1", default=False,
                        help="For L1 insertion calling, L1 specific function in clip-disc filtering")
    parser.add_argument("--gntp_feature",
                        action="store_true", dest="gntp_feature", default=False,
                        help="Collect genotyping features from bam")
    parser.add_argument("--gntp_classify",
                        action="store_true", dest="gntp_classify", default=False,
                        help="Train/predict genotpe classifier")
    parser.add_argument("--train_gntp",
                        action="store_true", dest="train_gntp", default=False,
                        help="Train the genotype classifer")
    parser.add_argument("--force",
                        action="store_true", dest="force", default=False,
                        help="Force to start from the very beginning")
    # parser.add_argument("--tumor",
    #                     action="store_true", dest="tumor", default=False,
    #                     help="Working on tumor samples")
####
    parser.add_argument("--bed",
                        action="store_true", dest="bed", default=False,
                        help="Input annotation in bed format")
    parser.add_argument("--flank", dest="flank", default=False,
                        help="flank regions")
    parser.add_argument("--gene",
                        action="store_true", dest="gene", default=False,
                        help="Check whether the insertion falls in genes")
    parser.add_argument("--user",
                        action="store_true", dest="user_specific", default=False,
                        help="User specific parameters, by default automatically calc the parameters")
    parser.add_argument("--single_sample",
                        action="store_true", dest="single_sample", default=False,
                        help="For single sample (like igv screenshot)")
    parser.add_argument("-i", "--input", dest="input", default="",
                        help="input file ", metavar="FILE")
    parser.add_argument("--input2", dest="input2", default="",
                        help="input file2 ", metavar="FILE")
    # parser.add_argument("-r", "--reference", dest="reference",
    #                     help="The list of files with repeat copies with flanking regions", metavar="FILE")
    parser.add_argument("--Alu-reference", dest="Alu_reference",
                        help="The file with Alu copies with flanking regions", metavar="FILE")
    parser.add_argument("--L1-reference", dest="L1_reference",
                        help="The file with L1 copies with flanking regions", metavar="FILE")
    parser.add_argument("--SVA-reference", dest="SVA_reference",
                        help="The file with SVA copies with flanking regions", metavar="FILE")
    # parser.add_argument("-a", "--annotation", dest="annotation",
    #                     help="The list of files with RepeatMasker annotations", metavar="FILE")
    parser.add_argument("--Alu-annotation", dest="Alu_annotation",
                        help="The file with Alu RepeatMasker annotations", metavar="FILE")
    parser.add_argument("--L1-annotation", dest="L1_annotation",
                        help="The file with L1 RepeatMasker annotations", metavar="FILE")
    parser.add_argument("--SVA-annotation", dest="SVA_annotation",
                        help="The file with SVA RepeatMasker annotations", metavar="FILE")
    # parser.add_argument("-c", "--copies", dest="copies",
    #                   help="Repeat copies ", metavar="FILE")
    parser.add_argument("-b", "--bam", dest="bam",
                        help="Input bam file", metavar="FILE")
    parser.add_argument("-d", "--barcode_bam", dest="barcode_bam",
                        help="Input barcode indexed bam file", metavar="FILE")
    parser.add_argument("-o", "--output", dest="output",
                        help="The output file", metavar="FILE")
    # YW 2021/05/10 added this
    # YW TO DO: add_subparsers to make this required for -D option
    parser.add_argument("-m", "--final_matrix", dest="final_matrix",
                        help="The final feature matrix", metavar="FILE")
    parser.add_argument("-p", "--path", dest="wfolder", type=str, default="./",
                        help="Working folder")
    parser.add_argument("--cp", dest="cwfolder", type=str,
                        help="Working folder for shared clipped reads")
    parser.add_argument("-n", "--cores", dest="cores", type=int, default=1,
                        help="number of cores")
    parser.add_argument("-e", "--extend", dest="extend", type=int, default=0,
                        help="extend length")
    parser.add_argument("--fflank", dest="fflank",
                        help="flank region file", metavar="FILE")
    parser.add_argument("--flklen", dest="flklen", type=int,
                        help="flank region file")
    # parser.add_argument("--purity", dest="purity", type=float, default=0.45,#by default tumor purity set to 45%
    #                     help="Tumor purity")
    parser.add_argument("--ref", dest="ref",
                        help="genome reference", metavar="FILE")
    # parser.add_argument("--cns", dest="cns",
    #                     help="list of files with repeat consensus", metavar="FILE") # YW need to change
    parser.add_argument("--Alu-cns", dest="Alu_cns",
                        help="Alu consensus file", metavar="FILE")
    parser.add_argument("--L1-cns", dest="L1_cns",
                        help="L1 consensus file", metavar="FILE")
    parser.add_argument("--SVA-cns", dest="SVA_cns",
                        help="SVA consensus file", metavar="FILE")
    parser.add_argument("--sc", dest="siteclip", type=int, default=2,
                        help="cutoff of minimum # of clipped reads at the exact position, use larger value for 10X")
    parser.add_argument("--lc", dest="lclip", type=int, default=3,
                        help="cutoff of minimum # of left clipped reads")
    parser.add_argument("--rc", dest="rclip", type=int, default=3,
                        help="cutoff of minimum # of right clipped reads")
    parser.add_argument("--cr", dest="cliprep", type=int, default=1,
                        help="cutoff of minimum # of clipped parts fall in repeats")
    parser.add_argument("--cc", dest="clipcns", type=int, default=1,
                        help="cutoff of minimum # of clipped parts fall in repeat cns")
    parser.add_argument("--nd", dest="ndisc", type=int, default=5,
                        help="cutoff of minimum # of discordant pair")
    # YW 2021/04/30 added this
    parser.add_argument("--ncns", dest="ncns", type=int, default=1,
                        help="cutoff of minimum # of clip + disc mapping to TE consensus")
    parser.add_argument("--cov", dest="cov", type=float, default=30.0,
                        help="approximate read depth")
    parser.add_argument("--iniclip", dest="iniclip", type=int, default=2,
                        help="initial minimum clip cutoff")
    parser.add_argument("--rmsk_extnd", dest="rmsk_extnd", type=int, default=100,
                        help="Length of the left extended region when loading the repeatmasker output")
    parser.add_argument("--rtype", dest="rep_type", type=int, default=1,
                        help="type of repeats: 1-L1, 2-Alu, 4-SVA, 8-HERV, 16-MIT, 32-MSTA")
    parser.add_argument("--blacklist", dest="blacklist", default="null",
                        help="Reference panel database for filtering, or a blacklist region", metavar="FILE")
    parser.add_argument("--model", dest="model", default="null",
                        help="Already trained model (.pkl file) for genotype classification", metavar="FILE")
    
    # YW 2021/05/26 added to parallelize cns remapping
    parser.add_argument("--c_realn_partition", dest="c_realn_partition", type=str, default="short", # medium for highcov
                        help="slurm partition for running realignment of clipped reads to repeat cns")
    parser.add_argument("--c_realn_time", dest="c_realn_time", type=str, default="0-08:00", #"1-00:00" for highcov
                        help="runtime for running realignment of clipped reads to repeat cns")
    parser.add_argument("--c_realn_mem", dest="c_realn_mem", type=int, default=50,
                        help="run memory (in GB) for running realignment of clipped reads to repeat cns")
    parser.add_argument("--d_realn_partition", dest="d_realn_partition", type=str, default="short",
                        help="slurm partition for running realignment of disc reads to repeat cns")
    parser.add_argument("--d_realn_time", dest="d_realn_time", type=str, default="0-00:10",
                        help="runtime for running realignment of disc reads to repeat cns")
    parser.add_argument("--d_realn_mem", dest="d_realn_mem", type=int, default=1,
                        help="run memory (in GB) for running realignment of disc reads to repeat cns")
    parser.add_argument("--check_interval", dest="check_interval", type=int, default=60,
                        help="interval (in sec) of checking whether the sbatch jobs of cns alignment of other repeat type have finished after finishing the main cns realignment job")
    # YW added 2023/03/29
    parser.add_argument("--email_user", dest="email_user", type=str, required=True,
                        help="specify the email address to notify SLURM job status")
    parser.add_argument("--sample_id", dest="sample_id", type=str, required=True,
                        help="specify sample id in job names for cns remapping and feature extraction, make it easy for job monitoring when running multiple samples")
    
    # YW 2021/07/27 added to enable control bam file coordinate lifting
    parser.add_argument("--ctrl", dest="ctrl",
                        action="store_true", default=False,
                        help="indicate running of aDNA ctrl bam input, skip feature extraction in -D step")
    parser.add_argument("--ref_bed", dest="ref_bed", default="/n/data1/bch/genetics/lee/elain/xTEA_benchmarking/alt_reference/gold_std_hg38_v4_int.temp",
                        help="The file with gold std set coordinate and insertion size, required only when --ctrl", metavar="FILE")
    parser.add_argument("--error_margin", dest="error_margin", type=int, default=15,
                        help="error margin to lift coordinate, required only when --ctrl")
    # YW 2021/07/30 added to enable input of ctrl bed file for coordinate subtraction
    parser.add_argument("--ctrl_bed", dest="ctrl_bed", default=None,
                        help="TEI coordinate from the control bam file, ancient only, required only when --ctrl==False", metavar="FILE")
    
    # YW 2021/09/29 added to enable extraction of clip info from pre-defined loci
    parser.add_argument("-L", "--locus_clip",
                        action="store_true", dest="locus_clip", default=False,
                        help="Extract clip info from predefined loci")
    parser.add_argument("--locus_file", dest="locus_file",
                        help="The bed file with predefined loci of polymorphic or reference TEIs, required only when -L", metavar="FILE")
    
    # YW 2021/10/27 added to enable extraction of features without a matched reference control
    parser.add_argument("--train",
                        action="store_true", dest="train", default=False,
                        help="Turn on training mode, where there is a matched reference control")
    
    # YW 2021/12/10 added to enable changing time limit for feature extraction by chunk
    parser.add_argument("--feat_extract_time", dest="feat_extract_time", default="0-01:30",
                        help="runtime for running feature extraction on 5M TEI sites")
    parser.add_argument("--feat_extract_partition", dest="feat_extract_partition", default="short",
                        help="partition for running feature extraction on 5M TEI sites")
    
    # YW 2022/10/03 added to enable merging features from the x_genotype_feature module
    parser.add_argument("--genotype_ft_output", dest="genotype_ft_output", default=None,
                        help="output from the x_genotype_feature module, input to the x_genotype_classify module when predicting")
    
    args = parser.parse_args()
    return args
####
####
# YW 2020/08/01 github update: add the last 2 arguments
# def automatic_gnrt_parameters(sf_bam_list, sf_ref, s_working_folder, n_jobs, b_force=False, b_tumor=False, f_purity=0.45):
def automatic_gnrt_parameters(sf_bam_list, sf_ref, s_working_folder, n_jobs, b_force=False):
    ####1. collect the basic information
    search_win = 500
    x_basic_info = X_BasicInfo(s_working_folder, n_jobs, sf_ref)
    rcd=x_basic_info.get_cov_is_rlth(sf_bam_list, sf_ref, search_win, b_force)
    f_cov=rcd[0]
    rlth=rcd[1]
    mean_is=rcd[2]
    std_var=rcd[3]

    ####2. based on the coverage, set the parameters
    xpar=Parameters()
    # if b_tumor==True:
    #     f_cov=f_cov*f_purity
    par_rcd=xpar.get_par_by_cov(f_cov) #in format (iclip, idisc, i_clip-disc)
    print("Ave coverage is {0}: automatic parameters (clip, disc, clip-disc) with value ({1}, {2} ,{3})\n".format(f_cov, par_rcd[0], par_rcd[1], par_rcd[2]))
    return par_rcd, rcd

####
def automatic_set_molecule_cutoff_for_10X_bam(sf_bam, sf_ref, s_working_folder, n_jobs):
    x_basic_info = X_BasicInfo(s_working_folder, n_jobs, sf_ref)
    search_win = 500
    f_cov=x_basic_info.calc_cov_from_bam(sf_bam, sf_ref, search_win)
    xpar=XParameters()
    f_molecule_cutoff=xpar.get_barcode_cov_cutoff(f_cov)
    return f_molecule_cutoff


def automatic_gnrt_parameters_case_control(sf_bam_list, sf_ref, s_working_folder, n_jobs, b_force=False):
    ####1. collect the basic information
    search_win = 500
    x_basic_info = X_BasicInfo(s_working_folder, n_jobs, sf_ref)
    rcd=x_basic_info.get_cov_is_rlth(sf_bam_list, sf_ref, search_win, b_force)
    f_cov=rcd[0]
    rlth=rcd[1]
    mean_is=rcd[2]
    std_var=rcd[3]
    ####2. based on the coverage, set the parameters
    xpar=CaseControlFilterPars()
    par_rcd=xpar.get_par_by_cov(f_cov) #in format (iclip, idisc, i_clip-disc)
    print("Ave coverage is {0}: automatic parameters (clip, disc, clip-disc) with value ({1}, {2} ,{3})\n".format(f_cov, par_rcd[0], par_rcd[1], par_rcd[2]))
    return par_rcd, rcd


# YW 2021/09/29 added for locus_clip option
# potentially add process_chrm_name
def load_locus_file(locus_file):
    m_list = {}
    with open(locus_file) as fin_candidate_sites:
        for line in fin_candidate_sites:
            if line[0] == "#":
                continue
            fields = line.split()
            if len(fields)<2: # chrm pos TE (optional)
                print(fields, " does not have enough fields")
                continue
            elif len(fields) == 2:
                chrm = fields[0]
                pos = int(fields[1])
                if chrm not in m_list:
                    m_list[chrm] = [pos]
                else:
                    m_list[chrm].append(pos)
            else: # 2021/11/29 added this to deal with input with both start and end positions
                chrm = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                if chrm not in m_list:
                    m_list[chrm] = [(start,end)]
                else:
                    m_list[chrm].append((start,end))
    return m_list

# YW 2021/12/03 added for dividing the candidate_list_from_disc.txt into equally sized chunks for memory issues
def count_lines(f_in):
    return int(subprocess.check_output(['wc', '-l', f_in]).split()[0])
    
####
# def adjust_cutoff_tumor(ncutoff=-1, i_adjust=1):
#     if ncutoff-i_adjust>1:
#         ncutoff=ncutoff-i_adjust
#     return ncutoff

####
##main function
if __name__ == '__main__':
    args = parse_arguments()

    if args.mit:#if this to call mitochondrial insertion, then will not filter out chrM in "x_intermediate_sites.py"
        global_values.turn_on_mit()

    if args.dna:
        global_values.turn_off_rna_mediated()
    if args.cbs:
        global_values.turn_on_check_by_sample()
    # if args.sva:
    #     global_values.turn_on_sva()
    # YW 2020/08/16 added this
    # if args.l1:
    #     global_values.turn_on_l1()

    b_automatic=True
    if args.user_specific:
        b_automatic=False
    # YW 2020/08/01 github update the following 3
    # b_tumor=args.tumor #whether this is tumor sample
    # f_purity=args.purity#tumor purity, by default 0.45
    b_resume=args.resume#resume the running, which will skip the step if output file already exist

    # YW 2023/03/29 added to customize email user for status update
    global_values.set_email_user(args.email_user)
    # for easy job monitoring
    global_values.set_sample_id(args.sample_id)

    
    if args.clip:  ###take in the normal illumina reads (10x will be viewed as normal illumina)
        print("Working on \"clip\" step!")
        sf_bam_list = args.input
        s_working_folder = args.wfolder # YW 2021/03/18 make this not specific to TE (since all processes are shared)
        n_jobs = args.cores
        sf_rep_cns_Alu =args.Alu_cns
        sf_rep_cns_L1 =args.L1_cns
        sf_rep_cns_SVA =args.SVA_cns
        sf_rep_Alu = args.Alu_reference  ####repeat copies "-r"
        sf_rep_L1 = args.L1_reference
        sf_rep_SVA = args.SVA_reference
        sf_annotation_Alu = args.Alu_annotation
        sf_annotation_L1 = args.L1_annotation
        sf_annotation_SVA = args.SVA_annotation
        
        sf_out = args.output
        b_se = args.single  ##single end reads or not, default is not
        sf_ref=args.ref ###reference genome "-ref"
        b_force=args.force #force to run from the very beginning
        # YW 2020/08/01 github update b_mosaic
        # b_mosaic=False #this is for mosaic calling from normal tissue # YW 2021/03/18 set this to False
        #i_iniclip=args.iniclip#
        if b_force == True:
            global_values.set_force_clean()
        site_clip_cutoff=args.siteclip #this is the cutoff for the exact position, use larger value for 10X
        global_values.set_initial_min_clip_cutoff(site_clip_cutoff)

        # merge the list from different bams of the same individual
        # Here when do the filtering, nearby regions are already considered!
        cutoff_left_clip = args.lclip
        cutoff_right_clip = args.rclip
        cutoff_clip_mate_in_rep = args.cliprep
        cutoff_clip_mate_in_cns = args.clipcns
        
        # YW 2021/05/26 added to parallelize cns remapping
        global_values.set_c_realign_partition(args.c_realn_partition)
        global_values.set_c_realign_time(args.c_realn_time)
        global_values.set_c_realign_memory(args.c_realn_mem)
        global_values.set_check_interval(args.check_interval)
        
        # YW 2020/08/01 github update: if statement and b_resume
        if b_resume == True and os.path.isfile(sf_out)==True:
            print("{0} exists, skipping \"clip\" step".format(sf_out))
        else:
            if os.path.isfile(sf_out)==True:
                print("User doesn't specify skipping, although {0} exists. Rerun the \"clip\" step.".format(sf_out))
            if b_automatic==True:
                # YW 2020/08/01 github update, 2 more arguments b_tumor, f_purity --> 2021/05/19 removed both
                rcd, basic_rcd=automatic_gnrt_parameters(sf_bam_list, sf_ref, s_working_folder, n_jobs, b_force)
                cutoff_left_clip=rcd[0]
                cutoff_right_clip=rcd[0]
                cutoff_clip_mate_in_rep=rcd[2]
            # YW 2020/06/28 added more annotations for the print message
            print("Clip cutoff: lclip: {0}, rclip: {1}, clip_mate_in_rep: {2}, clip_mate_in_cns: {3} are used!!!".format(cutoff_left_clip, cutoff_right_clip, cutoff_clip_mate_in_rep, cutoff_clip_mate_in_cns))
            tem_locator = TE_Multi_Locator(sf_bam_list, s_working_folder, n_jobs, sf_ref)

            ####by default, if number of clipped reads is larger than this value, then discard
            max_cov_cutoff=int(15*args.cov) #by default, this value is 600
            wfolder_pub_clip = args.cwfolder #public clip folder
            
            # YW 2020/08/09 clarified the comment
            ##Hard code inside:
            # 1. call_TEI_candidate_sites_from_clip_reads_v2 --> run_cnt_clip_part_aligned_to_rep_by_chrm_sort_version
            # here if 3/4 of the seq is mapped, then consider it as aligned to rep.
            ##2. require >=2 clip reads (actually depending on user input), whose clipped part is aligned to repeat copies (depending on --cr)
        
            # YW 2020/08/01 added b_mosaic input (github update), and updated in the function too, but this hasn't been used
            # YW 2021/05/19 took out b_mosaic
            tem_locator.call_TEI_candidate_sites_from_multiple_alignmts(sf_annotation_Alu, sf_annotation_L1, sf_annotation_SVA,
                                                                        sf_rep_cns_Alu, sf_rep_cns_L1, sf_rep_cns_SVA,
                                                                        sf_rep_Alu, sf_rep_L1, sf_rep_SVA,
                                                                        b_se, cutoff_left_clip,
                                                                        cutoff_right_clip, cutoff_clip_mate_in_rep, cutoff_clip_mate_in_cns,
                                                                        wfolder_pub_clip, b_force, max_cov_cutoff, sf_out)
####
    elif args.discordant:  # this views all the alignments as normal illumina reads
        print("Working on \"disc\" step!")
        sf_bam_list = args.bam  ###read in a bam list file
        s_working_folder = args.wfolder
        n_jobs = args.cores
        sf_rep_cns_Alu =args.Alu_cns
        sf_rep_cns_L1 =args.L1_cns
        sf_rep_cns_SVA =args.SVA_cns
        sf_annotation_Alu = args.Alu_annotation
        sf_annotation_L1 = args.L1_annotation
        sf_annotation_SVA = args.SVA_annotation
        sf_candidate_list = args.input # YW 2021/04/23 retained map counts to cns
        sf_out = args.output
        sf_ref = args.ref  ###reference genome, some cram file require this file to open
        
        # YW 2021/05/26 added to parallelize cns remapping
        global_values.set_d_realign_partition(args.d_realn_partition)
        global_values.set_d_realign_time(args.d_realn_time)
        global_values.set_d_realign_memory(args.d_realn_mem)
        global_values.set_check_interval(args.check_interval)
        
        feature_matrix = args.final_matrix # YW 2021/05/21 added this to output the final feature matrix
        peak_window = global_values.PEAK_WINDOW_DEFAULT # YW 2020/07/03 note: max distance between two clipped positions for them to be considered as from one insertion/cluster
        # if args.postFmosaic or args.somatic:#for mosaic events
        #     peak_window = global_values.PEAK_WINDOW_MOS_SOM
        # YW 2020/08/04 modified github update: originally will rerun if b_resume==False, now added extra if/elif cases
        if b_resume == True and os.path.isfile(sf_out)==True:
            print("{0} exists, skipping \"disc\" step".format(sf_out))
        else:
            if os.path.isfile(sf_out)==True:
                print("User doesn't specify skipping, although {0} exists. Rerun the \"disc\" step.".format(sf_out))
            xfilter = XIntermediateSites()
            m_original_sites = xfilter.load_in_candidate_list(sf_candidate_list)
            # YW CHANGED THE FUNCTION NAME from call_peak_candidate_sites_with_std_derivation
            # m_sites_clip_peak = xfilter.call_peak_candidate_sites_calc_std_deviation(m_original_sites, peak_window)
            m_sites_clip_peak = xfilter.call_peak_candidate_sites(m_original_sites, peak_window) # YW 2021/09/24 no longer need clip_pos_std here, so switch to a shorter function
            m_original_sites.clear()  #release the memory
            # YW added the following message
            print("Finished merging nearby clipped sites!")
            b_force = True
            rcd, basic_rcd = automatic_gnrt_parameters(sf_bam_list, sf_ref, s_working_folder, n_jobs, b_force)
            # YW TO DO: change to manual input of read length, insert size, and std deviation of insert size
            rlth = basic_rcd[1]  # read length
            mean_is = basic_rcd[2]  # mean insert size
            std_var = basic_rcd[3]  # standard deviation
            max_is = int(mean_is + 3 * std_var) + int(rlth) # max insert size, for search of discordant reads on the left/right of the clipped position
            iextend = max_is # see above
            
            i_is = global_values.INSERT_SIZE  ###set the insert size a large value, YW 2020/07/21 to tell whether two reads mapped to the same chr form a discordant read pair
            f_dev = std_var # --> what is this for?

            # this is the cutoff for  "left discordant" and "right discordant"
            # Either of them is larger than this cutoff, the site will be reported
            n_disc_cutoff = args.ndisc
            n_cns_cutoff = args.ncns # YW 2021/04/30 added this
            if b_automatic==True:
                n_disc_cutoff=rcd[1]

            print("Discordant cutoff: {0} is used!!!".format(n_disc_cutoff))

            # sf_tmp = s_working_folder + "disc_tmp.list" # YW 2021/04/30 commented this out because we don't need it
            # YW 2020/08/03 github update: added the following line
            sf_raw_disc=sf_out + global_values.RAW_DISC_TMP_SUFFIX #save the left and right raw disc for each site
            tem_locator = TE_Multi_Locator(sf_bam_list, s_working_folder, n_jobs, sf_ref)
            # YW 2020/08/03 github update: added 2 arguments sf_raw_disc, b_tumor --> 2021/05/19 removed b_tumor
            # YW 2021/04/30 removed sf_tmp from the input (no longer need this file output)
            tem_locator.filter_candidate_sites_by_discordant_pairs_multi_alignmts(m_sites_clip_peak, iextend, i_is, f_dev,
                                                                                  n_disc_cutoff,
                                                                                  sf_rep_cns_Alu, sf_rep_cns_L1, sf_rep_cns_SVA,
                                                                                  sf_annotation_Alu, sf_annotation_L1,
                                                                                  sf_annotation_SVA,
                                                                                  sf_raw_disc)
            m_sites_clip_peak.clear() # YW 2021/09/24 added to save memory
            
            # YW 2020/07/20 modified merge_clip_disc function to make sure locations without discordant read support will go through
            # xfilter.merge_clip_disc(sf_tmp, sf_candidate_list, sf_out)
            # YW 2021/04/21 wrote the function below to merge features from clip and disc
            xfilter.merge_clip_disc_new(sf_candidate_list, sf_raw_disc, sf_out + ".tmp", sf_out + ".clip_std_pos", n_cns_cutoff) # YW 2021/04/29 added the default cutoff of read count mapping to repeat cns
            # 2021/12/08 add clip_pos_std info to the sf_candidate_list output using bedtools intersect
            xfilter.add_clip_pos_std(sf_out + ".tmp", sf_out + ".clip_std_pos", sf_out)
        
        cmd_runner = CMD_RUNNER()
        b_train=args.train
        if args.train: # YW 2021/10/27 added this if statement
            if args.ctrl == False:
                if b_resume and os.path.isfile(feature_matrix):
                    if os.path.getsize(feature_matrix)>0:
                        print(f"{feature_matrix} exists, skipping \"feature extraction\" step!")
                    else:
                        coor_lift = Coor_Lift(sf_out, sf_out + ".sorted", None, args.error_margin)
                        coor_lift.sort_subtract_overlap(args.ctrl_bed)
                        feat_folder = s_working_folder + "features/"
                        feat_mat = Feature_Matrix(sf_out + ".sorted", feature_matrix+".tmp", feat_folder, sf_bam_list, sf_ref, n_jobs, b_train)
                        feat_mat.run_feature_extraction()
                else:
                    if b_resume and os.path.isfile(sf_out + ".sorted"):
                        if os.path.getsize(sf_out + ".sorted")>0:
                            print(f"{sf_out}.sorted exists. Proceed to generating the feature matrix...")
                        else:
                            coor_lift = Coor_Lift(sf_out, sf_out + ".sorted", None, args.error_margin)
                            coor_lift.sort_subtract_overlap(args.ctrl_bed)
                    else:
                        # YW 2021/07/30 wrote the function below to subtract TEI coordinates overlapping with ctrl (ancient only!!!)
                        coor_lift = Coor_Lift(sf_out, sf_out + ".sorted", None, args.error_margin)
                        coor_lift.sort_subtract_overlap(args.ctrl_bed)
                    # YW 2021/05/10 wrote the function below to extract extra features
                    feat_folder = s_working_folder + "features/"
                    feat_mat = Feature_Matrix(sf_out + ".sorted", feature_matrix+".tmp", feat_folder, sf_bam_list, sf_ref, n_jobs, b_train)
                    feat_mat.run_feature_extraction()
                cmd_runner.run_cmd_to_file(f"sort -V -k 1,2 {feature_matrix}.tmp", feature_matrix)
                os.remove(feature_matrix + ".tmp")
            else: # YW 2021/07/27 perform coordinate lifting to gold std set deleted ref
                if b_resume and os.path.isfile(sf_out + ".lifted"):
                    if os.path.getsize(sf_out + ".lifted")>0:
                        print(f"{sf_out}.lifted exists, skipping \"coordinate lifting\" step for ctrl!")
                    else:
                        coor_lift = Coor_Lift(sf_out, sf_out + ".lifted", args.ref_bed, args.error_margin)
                        coor_lift.run_coor_lift("ctrl")
                else:
                    coor_lift = Coor_Lift(sf_out, sf_out + ".lifted", args.ref_bed, args.error_margin)
                    coor_lift.run_coor_lift("ctrl")
        else: # YW 2021/10/27 added the following
            global_values.set_feat_extract_time(args.feat_extract_time)
            global_values.set_feat_extract_partition(args.feat_extract_partition)
            if b_resume and os.path.isfile(feature_matrix):
                if os.path.getsize(feature_matrix) > 0:
                    sys.exit(f"{feature_matrix} exists. Exiting...")
            
            # YW 2021/12/03 changed the following to break the sf_out into chunks (for potential parallelization)
            line_cnt = count_lines(sf_out)
            idx_lst = []
            if line_cnt < global_values.CHUNK_SIZE: # no need to generate chunks for small input
                # cmd_runner.run_cmd_to_file(f"sort -V -k 1,2 {sf_out}", sf_out + ".sorted")
                feat_folder = s_working_folder + "features/"
                feat_mat = Feature_Matrix(sf_out, feature_matrix, feat_folder, sf_bam_list, sf_ref, n_jobs, b_train) 
                feat_mat.run_feature_extraction()
            else:
                feat_prefix = sf_out + "_"
                feat_suffix = ".tmp"
                split_disc = f"shuf {sf_out} | split -l {global_values.CHUNK_SIZE} -d --additional-suffix=.tmp - {feat_prefix}"
                cmd_runner.run_cmd_small_output(split_disc)
                list_sf_out = glob(f"{feat_prefix}*{feat_suffix}")
                #######################################
                # YW 2021/12/07 to launch multiple sbatch jobs for all chunks of feature extraction
                # currently, this hasn't been tested for when there is more than 1 bam file
                feat_folder = s_working_folder + "features/"
                feat_parallel = FEAT_PARALLEL(s_working_folder, feat_folder, sf_bam_list, sf_ref, n_jobs, b_train)
                sbatch_scr_list = []
                for tmp_file in list_sf_out[1:]:
                    f_name = os.path.basename(tmp_file)
                    idx = f_name.lstrip(f"{os.path.basename(sf_out)}_").rstrip(feat_suffix)
                    idx_lst.append(idx)
                    py_script, done_f, fail_f = feat_parallel.gnrt_feat_py_scripts(tmp_file, feature_matrix, idx)
                    sbatch_scr_list.append((py_script, done_f, fail_f))
                
                feat_parallel.run_sbatch_scripts(sbatch_scr_list, idx_lst)
                
                tmp_file = list_sf_out[0]
                f_name = os.path.basename(tmp_file)
                idx = f_name.lstrip(f"{os.path.basename(sf_out)}_").rstrip(feat_suffix)
                idx_lst.append(idx)
                feat_mat = Feature_Matrix(tmp_file, feature_matrix+idx, feat_folder, sf_bam_list, sf_ref, n_jobs, b_train, idx)
                feat_mat.run_feature_extraction()
                print(f"{idx} feature extraction has finished")
                
                num_finished = 1 + len(glob(s_working_folder + "*_feat_extract.done"))
                num_failed = len(glob(s_working_folder + "*_feat_extract.fail"))
                if num_failed > 0:
                    sys.exit("Exit... One of the feature extraction jobs failed!!!")
                while num_finished < len(idx_lst):
                    sleep(global_values.CHECK_INTERVAL)
                    num_finished = 1 + len(glob(s_working_folder + "*_feat_extract.done"))
                    num_failed = len(glob(s_working_folder + "*_feat_extract.fail"))
                    if num_failed > 0:
                        sys.exit("Exit... One of the feature extraction jobs failed!!!")
                ########################################
                # concatenate and sort the feature matrix output
                cmd_runner.run_cmd_to_file(f"cat {feature_matrix}* | sort -V -k 1,2", feature_matrix)
                # remove intermediate files if sf_out has been split
                if idx_lst:
                    for idx in idx_lst:
                        os.remove(f"{sf_out}_{idx}.tmp")
                        os.remove(feature_matrix + idx)
                # remove the .done files
                sbatch_stat_list = glob(s_working_folder + "*_feat_extract.*")
                if sbatch_stat_list:
                    for f in sbatch_stat_list:
                        os.remove(f)
    
    # 2021/09/29 NEW OPTIONS FOR PRE-DEFINED LOCI
    elif args.locus_clip:
        print("Working on \"collecting clip info based on loci\" step!")
        locus_dict = load_locus_file(args.locus_file)
        sf_bam_list = args.input
        # UPDATE THE FOLLOWING in cmd!!!
        s_working_folder = args.wfolder # YW 2021/03/18 make this not specific to TE (since all processes are shared)
        n_jobs = args.cores
        sf_rep_cns_Alu =args.Alu_cns
        sf_rep_cns_L1 =args.L1_cns
        sf_rep_cns_SVA =args.SVA_cns
        sf_rep_Alu = args.Alu_reference  ####repeat copies "-r"
        sf_rep_L1 = args.L1_reference
        sf_rep_SVA = args.SVA_reference
        sf_annotation_Alu = args.Alu_annotation
        sf_annotation_L1 = args.L1_annotation
        sf_annotation_SVA = args.SVA_annotation
        
        sf_out = args.output
        b_se = args.single  ##single end reads or not, default is not
        sf_ref=args.ref ###reference genome "-ref"
        b_force=args.force #force to run from the very beginning
        # YW 2020/08/01 github update b_mosaic
        # b_mosaic=False #this is for mosaic calling from normal tissue # YW 2021/03/18 set this to False
        #i_iniclip=args.iniclip#
        if b_force == True:
            global_values.set_force_clean()
        site_clip_cutoff=args.siteclip #this is the cutoff for the exact position, use larger value for 10X
        global_values.set_initial_min_clip_cutoff(site_clip_cutoff)
        
        # YW 2021/09/29 should we keep these?
        # merge the list from different bams of the same individual
        # Here when do the filtering, nearby regions are already considered!
        cutoff_left_clip = args.lclip
        cutoff_right_clip = args.rclip
        cutoff_clip_mate_in_rep = args.cliprep
        cutoff_clip_mate_in_cns = args.clipcns
        
        # YW 2021/05/26 added to parallelize cns remapping
        global_values.set_c_realign_partition(args.c_realn_partition)
        global_values.set_c_realign_time(args.c_realn_time)
        global_values.set_c_realign_memory(args.c_realn_mem)
        global_values.set_check_interval(args.check_interval)
        
        # YW 2020/08/01 github update: if statement and b_resume
        if b_resume == True and os.path.isfile(sf_out)==True:
            print("{0} exists, skipping \"clip\" step".format(sf_out))
        else:
            if os.path.isfile(sf_out)==True:
                print("User doesn't specify skipping, although {0} exists. Rerun the \"clip\" step.".format(sf_out))
            if b_automatic==True:
                # YW 2020/08/01 github update, 2 more arguments b_tumor, f_purity --> 2021/05/19 removed both
                rcd, basic_rcd=automatic_gnrt_parameters(sf_bam_list, sf_ref, s_working_folder, n_jobs, b_force)
                cutoff_left_clip=rcd[0]
                cutoff_right_clip=rcd[0]
                cutoff_clip_mate_in_rep=rcd[2]
            # YW 2020/06/28 added more annotations for the print message
            print("Clip cutoff: lclip: {0}, rclip: {1}, clip_mate_in_rep: {2}, clip_mate_in_cns: {3} are used!!!".format(cutoff_left_clip, cutoff_right_clip, cutoff_clip_mate_in_rep, cutoff_clip_mate_in_cns))
            tem_locator = TE_Multi_Locator(sf_bam_list, s_working_folder, n_jobs, sf_ref)

            ####by default, if number of clipped reads is larger than this value, then discard
            max_cov_cutoff=int(15*args.cov) #by default, this value is 600
            wfolder_pub_clip = args.cwfolder #public clip folder (for predefined loci)
            
            # YW 2020/08/09 clarified the comment
            ##Hard code inside:
            # 1. call_TEI_candidate_sites_from_clip_reads_v2 --> run_cnt_clip_part_aligned_to_rep_by_chrm_sort_version
            # here if 3/4 of the seq is mapped, then consider it as aligned to rep.
            ##2. require >=2 clip reads (actually depending on user input), whose clipped part is aligned to repeat copies (depending on --cr)
        
            # YW 2020/08/01 added b_mosaic input (github update), and updated in the function too, but this hasn't been used
            # YW 2021/05/19 took out b_mosaic
            tem_locator.collect_clip_info_from_TEI_candidate_sites(locus_dict,
                                                                   sf_annotation_Alu, sf_annotation_L1, sf_annotation_SVA,
                                                                   sf_rep_cns_Alu, sf_rep_cns_L1, sf_rep_cns_SVA,
                                                                   sf_rep_Alu, sf_rep_L1, sf_rep_SVA,
                                                                   b_se, cutoff_left_clip,
                                                                   cutoff_right_clip, cutoff_clip_mate_in_rep, cutoff_clip_mate_in_cns,
                                                                   wfolder_pub_clip, b_force, max_cov_cutoff, sf_out)
####
####
    ####
    elif args.filter_csn:  #filter out the FP by the pattern in the consensus repeat
        print("Working on \"clip-disc-filtering\" step!")
        sf_bam_list = args.bam  ###read in a bam list file
        s_working_folder = args.wfolder
        if s_working_folder[-1] != "/":
            s_working_folder += "/"
        print("Current working folder is: {0}\n".format(s_working_folder))
        n_jobs = args.cores
        sf_ref = args.ref  ###reference genome, some cram file require this file to open

        sf_candidate_list = args.input
        sf_raw_disc=args.input2#this is the raw disc file
        # YW 2020/07/05 changed the following from 400 to 200 (shorter read length in ancient samples)
        iextnd = global_values.F_CNS_EXTEND ###for each site, re-collect reads in range [-iextnd, iextnd], this around ins +- 3*deviation
        bin_size = global_values.BIN_SIZE  # block size for parallelization
        sf_rep_cns_Alu =args.Alu_cns
        sf_rep_cns_L1 =args.L1_cns
        sf_rep_cns_SVA =args.SVA_cns
        bmapped_cutoff = global_values.MIN_CLIP_MAPPED_RATIO # minimal ratio of aligned bases in clipped part to be qualified (in both alignment of clipped part in clipped and in discordant to repeat cns)
        i_concord_dist = 550  # this should be the mean_is+3*is_std_deviation, used to cluster disc reads on the consensus YW: --> global_values? (too large?) YW 2020/07/19 clarified this comment
        f_concord_ratio = global_values.DISC_CONCORD_RATIO # YW 2020/07/21: to check whether discordant reads are clustered on cns, currently disabled
        sf_output = args.output
        sf_flank=args.fflank
        i_flank_lenth = args.flklen
        
        # YW 2020/08/04 modified github update: originally will rerun if b_resume==False, now added extra if/elif cases
        if b_resume == True and os.path.isfile(sf_output)==True:
            print("{0} exists, skipping \"clip-disc-filtering\" step.\n".format(sf_output))
        else:
            if os.path.isfile(sf_output)==True:
                print("User doesn't specify skipping, although {0} exists. Rerun the \"clip-disc-filtering\" step.\n".format(sf_output))
            b_force=True
            rcd, basic_rcd = automatic_gnrt_parameters(sf_bam_list, sf_ref, s_working_folder, n_jobs, b_force)
            ave_cov = basic_rcd[0]  # ave coverage
            rlth = basic_rcd[1]  # read length
            mean_is = basic_rcd[2]  # mean insert size
            std_var = basic_rcd[3]  # standard deviation
            print("Mean insert size is: {0}\n".format(mean_is))
            print("Standard deviation is: {0}\n".format(std_var))
            
            max_is = int(mean_is + 3 * std_var)
            if iextnd < max_is: #correct the bias
                iextnd = max_is
            if i_concord_dist < max_is: #correct the bias
                i_concord_dist = max_is
            global_values.set_read_length(rlth)
            global_values.set_insert_size(max_is)
            global_values.set_average_cov(ave_cov)
            print("Read length is: {0}\n".format(rlth))
            print("Maximum insert size is: {0}\n".format(max_is))
            print("Average coverage is: {0}\n".format(ave_cov))

            n_clip_cutoff = args.cliprep #this is the sum of left and right clipped reads
            n_disc_cutoff = args.ndisc  #each sample should have at least this number of discordant reads
            if b_automatic==True:
                n_clip_cutoff=rcd[0]
                n_disc_cutoff=rcd[1]
            print("Filter (on cns) cutoff: number of clipped: {0} and number of discordant reads: {1} are used!!!\n".format(n_clip_cutoff, n_disc_cutoff))

            x_cd_filter = XClipDiscFilter(sf_bam_list, s_working_folder, n_jobs, sf_ref)
            x_cd_filter.call_MEIs_consensus(sf_candidate_list, iextnd, bin_size,
                                            sf_rep_cns_Alu, sf_rep_cns_L1, sf_rep_cns_SVA,
                                            sf_flank, i_flank_lenth,
                                            bmapped_cutoff, i_concord_dist, f_concord_ratio, n_clip_cutoff, n_disc_cutoff,
                                            sf_output)
    elif args.gntp_feature:#generate the genotype features
        sf_bam_list = args.bam
        sf_ref = args.ref
        sf_candidate_list = args.input
        n_jobs = args.cores
        s_working_folder = args.wfolder
        sf_output = args.output

        x_gntper = XGenotyper(sf_ref, s_working_folder, n_jobs)
        extnd = 450
        x_gntper.call_genotype(sf_bam_list, sf_candidate_list, extnd, sf_output)

    elif args.gntp_classify:
        b_train = args.train_gntp
        sf_model=args.model
        if b_train==True:#train a new model
            print("training a new model")
            # sf_00_list = "/n/data1/hms/dbmi/park/simon_chu/projects/XTEA/genotyping/training_set_SSC/Genotyping/rslt_list/all_00.list"
            # sf_01_list = "/n/data1/hms/dbmi/park/simon_chu/projects/XTEA/genotyping/training_set_SSC/Genotyping/rslt_list/all_01.list"
            # sf_11_list = "/n/data1/hms/dbmi/park/simon_chu/projects/XTEA/genotyping/training_set_SSC/Genotyping/rslt_list/all_11.list"
            # sf_arff = "/n/data1/hms/dbmi/park/simon_chu/projects/XTEA/genotyping/training_set_SSC/Genotyping/merged_all_0_1_2.arff"
            sf_01_list = "/n/data1/bch/genetics/lee/elain/hapROH/xTea_genotyper/train_files/all_01.list"
            sf_11_list = "/n/data1/bch/genetics/lee/elain/hapROH/xTea_genotyper/train_files/all_11.list"
            sf_arff = args.input
            gc = GntpClassifier_DF21()
            b_balance=False
            gc.gnrt_training_arff_from_xTEA_output(sf_01_list, sf_11_list, sf_arff, b_balance)
            #pkl_filename = "./genotyping/trained_model_ssc_py2_random_forest_two_category.pkl"
            gc.train_model(sf_arff, sf_model)
        else:#predict the genotype
            #sf_model = "./genotyping/trained_model_ssc_py2_random_forest_two_category.pkl"
            sf_xTEA = args.input #input raw results before calling genotype
            sf_new = args.output
            sf_genotype_ft = args.genotype_ft_output
            gc = GntpClassifier_DF21()
            # pkl_model = gc.load_model_from_file(sf_model)
            # sf_arff = sf_xTEA + ".arff"
            gc.predict_for_site(sf_model, sf_xTEA, sf_genotype_ft, sf_new)

