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
import global_values
from x_TEI_locator import *
# from x_TEI_source_tracer import *
from x_local_assembly import *
from x_intermediate_sites import *
from x_reference import *
from x_clip_disc_filter import *
from x_somatic_calling import *
# from x_analysis import *
import argparse
from x_reads_collection import *
from x_mutation import *
from x_sv import *
from x_gene_annotation import *
from x_genotype_feature import *
from x_basic_info import *
from x_parameter import *
from x_post_filter import *
from x_mosaic_calling import *
from x_joint_calling import *
from x_igv import *
from x_gvcf import * # new
from x_genotype_classify import * # new
from x_orphan_transduction import * # new

####
##parse the options
def parse_arguments():
    parser = argparse.ArgumentParser("Main script of xTea_ML")
    parser.add_argument("-P", "--preprocess",
                        action="store_true", dest="preprocess", default=False,
                        help="Preprocessing stpes")
    parser.add_argument("-Q","--collectclip",
                        action="store_true", dest="collect_clip", default=False,
                        help="Call clipped reads from alignment")
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
    #transduction
    parser.add_argument("--transduction",
                        action="store_true", dest="transduction", default=False,
                        help="Call transduction for sites")
    parser.add_argument("--sibling",
                        action="store_true", dest="sibling", default=False,
                        help="Call sibling transduction for sites")
    parser.add_argument("--spectrum",
                        action="store_true", dest="spectrum", default=False,
                        help="Spectrum analysis by tumor type")
    parser.add_argument("-B", "--barcode",
                        action="store_true", dest="barcode", default=False,
                        help="Indicate the input is 10X bam")
    parser.add_argument("-E", "--collect",
                        action="store_true", dest="collect", default=False,
                        help="Collect reads for candidate sites")
    parser.add_argument("-I", "--mutation",
                        action="store_true", dest="mutation", default=False,
                        help="Call internal mutation")
    parser.add_argument("-U", "--collect_Illumina",
                        action="store_true", dest="collect_illumina", default=False,
                        help="Collect reads for candidate sites from normal illumina alignment")
    parser.add_argument("-F", "--filter_asm",
                        action="store_true", dest="filter_asm", default=False,
                        help="Filter out candidate sites from assembly")
    parser.add_argument("-G", "--contig_realign",
                        action="store_true", dest="contig_realign", default=False,
                        help="Filter out candidate sites from assembly")
    parser.add_argument("-T", "--trace",
                        action="store_true", dest="trace", default=False,
                        help="Trace the sources of TEIs")
    parser.add_argument("-A", "--assembly",
                        action="store_true", dest="assembly", default=False,
                        help="Do local assembly for collected reads")
    parser.add_argument("-L", "--local",
                        action="store_true", dest="local", default=False,
                        help="Assemble the TEIs on local machine")
    parser.add_argument("-M", "--map",
                        action="store_true", dest="map", default=False,
                        help="map flank regions to the assembled contigs")
    parser.add_argument("-V", "--visualization",
                        action="store_true", dest="visualization", default=False,
                        help="Show the heatmap figure of the selected regions")
    parser.add_argument("-K", "--withflank",
                        action="store_true", dest="withflank", default=False,
                        help="Keep the flank regions with the repeat copies")
    parser.add_argument("-J", "--joint",
                        action="store_true", dest="joint", default=False,
                        help="Joint calling")
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
    parser.add_argument("--postF",
                        action="store_true", dest="postF", default=False,
                        help="Post filtering module")
    parser.add_argument("--gntp_classify",
                        action="store_true", dest="gntp_classify", default=False,
                        help="Train/predict genotpe classifier")
    parser.add_argument("--train_gntp",
                        action="store_true", dest="train_gntp", default=False,
                        help="Train the genotype classifer")
    parser.add_argument("--postFmosaic",
                        action="store_true", dest="postFmosaic", default=False,
                        help="Post filtering module for mosaic events")
    parser.add_argument("--igv",
                        action="store_true", dest="igv", default=False,
                        help="Prepare screenshot command for given sites")
    parser.add_argument("--force",
                        action="store_true", dest="force", default=False,
                        help="Force to start from the very beginning")
    parser.add_argument("--case_control",
                        action="store_true", dest="case_control", default=False,
                        help="case control mode")
    parser.add_argument("--tumor",
                        action="store_true", dest="tumor", default=False,
                        help="Working on tumor samples")
    #convert to gVCF
    parser.add_argument("--gVCF",
                        action="store_true", dest="gVCF", default=False,
                        help="Generate the gVCF from xTEA raw output")
####
    parser.add_argument("--bed",
                        action="store_true", dest="bed", default=False,
                        help="Input annotation in bed format")
    parser.add_argument("--mosaic",
                        action="store_true", dest="mosaic", default=False,
                        help="Call mosaic events")
    parser.add_argument("--flk_map",
                        action="store_true", dest="flk_map", default=False,
                        help="Map flanks to contigs")
    parser.add_argument("--analysis",
                        action="store_true", dest="analysis", default=False,
                        help="Result analysis")
    parser.add_argument("--flank", dest="flank", default=False,
                        help="flank regions")
    parser.add_argument("--sv",
                        action="store_true", dest="sv", default=False,
                        help="Call promoted SVs")
    parser.add_argument("--gene",
                        action="store_true", dest="gene", default=False,
                        help="Check whether the insertion falls in genes")
    parser.add_argument("--somatic",
                        action="store_true", dest="somatic", default=False,
                        help="Only call somatic events from high coverage normal samples")
    parser.add_argument("--somatic_hc",
                        action="store_true", dest="somatic_hc", default=False,
                        help="Get high confident somatic events")
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
    parser.add_argument("-p", "--path", dest="wfolder", type="string", default="./",
                        help="Working folder")
    parser.add_argument("--cp", dest="cwfolder", type="string",
                        help="Working folder for shared clipped reads")
    parser.add_argument("-n", "--cores", dest="cores", type="int", default=1,
                        help="number of cores")
    parser.add_argument("-e", "--extend", dest="extend", type="int", default=0,
                        help="extend length")
    parser.add_argument("-u", "--dup", dest="duplication",
                        help="duplication files", metavar="FILE")
    parser.add_argument("--fflank", dest="fflank",
                        help="flank region file", metavar="FILE")
    parser.add_argument("--flklen", dest="flklen", type="int",
                        help="flank region file")
    parser.add_argument("--purity", dest="purity", type="float", default=0.45,#by default tumor purity set to 45%
                        help="Tumor purity")
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
    parser.add_argument("--sc", dest="siteclip", type="int", default=2,
                        help="cutoff of minimum # of clipped reads at the exact position, use larger value for 10X")
    parser.add_argument("--lc", dest="lclip", type="int", default=3,
                        help="cutoff of minimum # of left clipped reads")
    parser.add_argument("--rc", dest="rclip", type="int", default=3,
                        help="cutoff of minimum # of right clipped reads")
    parser.add_argument("--cr", dest="cliprep", type="int", default=1,
                        help="cutoff of minimum # of clipped parts fall in repeats")
    parser.add_argument("--nd", dest="ndisc", type="int", default=5,
                        help="cutoff of minimum # of discordant pair")
    parser.add_argument("--nb", dest="nbarcode", type="int", default=500,
                        help="cutoff of maximum # of molecure coverage")
    parser.add_argument("--teilen", dest="teilen", type="int",
                        help="minimum length of the insertion for future analysis")
    parser.add_argument("--cov", dest="cov", type="float", default=30.0,
                        help="approximate read depth")
    parser.add_argument("--iniclip", dest="iniclip", type="int", default=2,
                        help="initial minimum clip cutoff")
    parser.add_argument("--af1", dest="af1", type="float", default=0.005,
                        help="minimal allel fraction")
    parser.add_argument("--af2", dest="af2", type="float", default=0.45,
                        help="minimal allel fraction")
    parser.add_argument("--rmsk_extnd", dest="rmsk_extnd", type="int", default=100,
                        help="Length of the left extended region when loading the repeatmasker output")
    parser.add_argument("--rtype", dest="rep_type", type="int", default=1,
                        help="type of repeats: 1-L1, 2-Alu, 4-SVA, 8-HERV, 16-MIT, 32-MSTA")
    parser.add_argument("--blacklist", dest="blacklist", default="null",
                        help="Reference panel database for filtering, or a blacklist region", metavar="FILE")
    parser.add_argument("--model", dest="model", default="null",
                        help="Already trained model (.pkl file) for genotype classification", metavar="FILE")
    args = parser.parse_args()
    return args
####
####
# YW 2020/08/01 github update: add the last 2 arguments
def automatic_gnrt_parameters(sf_bam_list, sf_ref, s_working_folder, n_jobs, b_force=False, b_tumor=False, f_purity=0.45):
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
    if b_tumor==True:
        f_cov=f_cov*f_purity
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

####
def adjust_cutoff_tumor(ncutoff=-1, i_adjust=1):
    if ncutoff-i_adjust>1:
        ncutoff=ncutoff-i_adjust
    return ncutoff

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
    if args.sva:
        global_values.turn_on_sva()
    # YW 2020/08/16 added this
    if args.l1:
        global_values.turn_on_l1()

    b_automatic=True
    if args.user_specific:
        b_automatic=False
    # YW 2020/08/01 github update the following 3
    b_tumor=args.tumor #whether this is tumor sample
    f_purity=args.purity#tumor purity, by default 0.45
    b_resume=args.resume#resume the running, which will skip the step if output file already exist

    if args.preprocess:  #preprocess steps
        s_working_folder = args.wfolder
        sf_ref = args.reference
        sf_annotation = args.annotation
        sf_out_fa = args.output
        flank_lth = args.extend
        b_with_flank = args.withflank  # if not set then no flank region
        b_bed_fmt=args.bed
        x_annotation = XAnnotation(sf_annotation)
        i_lextnd=args.rmsk_extnd
        global_values.set_load_rmsk_left_extnd(i_lextnd)
        if b_bed_fmt==True:
            x_annotation.collect_seqs_of_TE_from_ref_bed_fmt(sf_ref, sf_out_fa, flank_lth)
        else:# this is for repeatmasker output
            b_with_chr = x_annotation.is_ref_chrm_with_chr(sf_ref)
            x_annotation.set_with_chr(b_with_chr)  # if chrm in reference has "chr", then True, otherwise False
            x_annotation.load_rmsk_annotation()
            # x_annotation.collect_flank_regions_of_TE_from_ref(sf_ref, flank_lth, sf_out_fa) #only get the flank regions
            x_annotation.collect_seqs_of_TE_from_ref(sf_ref, flank_lth, b_with_flank, sf_out_fa)
            x_annotation.bwa_index_TE_seqs(sf_out_fa)
    elif args.flank:#preprocess the flank regions steps
        s_working_folder = args.wfolder
        sf_ref = args.reference
        sf_annotation = args.annotation
        sf_out_fa = args.output
        flank_lth = args.extend
        x_annotation = XAnnotation(sf_annotation)
        b_with_chr = x_annotation.is_ref_chrm_with_chr(sf_ref)
        x_annotation.set_with_chr(b_with_chr)  # if chrm in reference has "chr", then True, otherwise False
        b_bed_fmt = args.bed
        if b_bed_fmt==True:
            print("load from bed")
            x_annotation.load_annotation_no_extnd_from_bed()
        else:
            print("load from rmsk")
            x_annotation.load_rmsk_annotation_no_extnd()
        x_annotation.collect_flank_regions_of_TE_from_ref(sf_ref, flank_lth, sf_out_fa)  # only get the flank regions
        x_annotation.bwa_index_TE_seqs(sf_out_fa)


####
    elif args.clip:  ###take in the normal illumina reads (10x will be viewed as normal illumina)
        print("Working on \"clip\" step!")
        sf_bam_list = args.input
        s_working_folder = args.wfolder
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
        b_mosaic=args.mosaic #this is for mosaic calling from normal tissue
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
        
        # YW 2020/08/01 github update: if statement and b_resume
        # YW 2020/08/04 modified github update: originally will rerun if b_resume==False, now added extra if/elif cases
        if b_resume == True and os.path.isfile(sf_out)==True:
            print("{0} exists, skipping \"clip\" step".format(sf_out))
        else:
            if os.path.isfile(sf_out)==True:
                print("User doesn't specify skipping, although {0} exists. Rerun the \"clip\" step.".format(sf_out))
            if b_automatic==True:
                # YW 2020/08/01 github update, 2 more arguments b_tumor, f_purity
                rcd, basic_rcd=automatic_gnrt_parameters(sf_bam_list, sf_ref, s_working_folder, n_jobs,
                                                         b_force, b_tumor, f_purity)
                cutoff_left_clip=rcd[0]
                cutoff_right_clip=rcd[0]
                # if b_tumor==True:
                #     cutoff_left_clip=adjust_cutoff_tumor(cutoff_left_clip)
                #     cutoff_right_clip=adjust_cutoff_tumor(cutoff_right_clip)
                cutoff_clip_mate_in_rep=rcd[2]
            # YW 2020/06/28 added more annotations for the print message
            print("Clip cutoff: lclip: {0}, rclip: {1}, clip_mate_in_rep: {2} are used!!!".format(cutoff_left_clip, cutoff_right_clip, cutoff_clip_mate_in_rep))
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
            tem_locator.call_TEI_candidate_sites_from_multiple_alignmts(sf_annotation_Alu, sf_annotation_L1, sf_annotation_SVA,
                                                                        sf_rep_cns_Alu, sf_rep_cns_L1, sf_rep_cns_SVA,
                                                                        sf_rep_Alu, sf_rep_L1, sf_rep_SVA,
                                                                        b_se, cutoff_left_clip,
                                                                        cutoff_right_clip, cutoff_clip_mate_in_rep, b_mosaic,
                                                                        wfolder_pub_clip, b_force, max_cov_cutoff, sf_out)
####
    elif args.discordant:  # this views all the alignments as normal illumina reads
        print("Working on \"disc\" step!")
        sf_bam_list = args.bam  ###read in a bam list file
        s_working_folder = args.wfolder
        n_jobs = args.cores
        sf_annotation_Alu = args.Alu_annotation
        sf_annotation_L1 = args.L1_annotation
        sf_annotation_SVA = args.SVA_annotation
        sf_candidate_list = args.input
        sf_out = args.output
        sf_ref = args.ref  ###reference genome, some cram file require this file to open
        peak_window = global_values.PEAK_WINDOW_DEFAULT # YW 2020/07/03 note: max distance between two clipped positions for them to be considered as from one insertion/cluster
        if args.postFmosaic or args.somatic:#for mosaic events
            peak_window = global_values.PEAK_WINDOW_MOS_SOM
        # YW 2020/08/04 modified github update: originally will rerun if b_resume==False, now added extra if/elif cases
        if b_resume == True and os.path.isfile(sf_out)==True:
            print("{0} exists, skipping \"disc\" step".format(sf_out))
        else:
            if os.path.isfile(sf_out)==True:
                print("User doesn't specify skipping, although {0} exists. Rerun the \"disc\" step.".format(sf_out))
            xfilter = XIntermediateSites()
            m_original_sites = xfilter.load_in_candidate_list(sf_candidate_list)
            sf_peak_sites = s_working_folder + "clip_peak_candidate.list"
            #m_sites_clip_peak = xfilter.call_peak_candidate_sites(m_original_sites, PEAK_WINDOW)  # get the peak sites
                # get the peak sites
            # YW CHANGED THE FUNCTION NAME from call_peak_candidate_sites_with_std_derivation
            m_sites_clip_peak = xfilter.call_peak_candidate_sites_calc_std_deviation(m_original_sites, peak_window)
            xfilter.output_candidate_sites(m_sites_clip_peak, sf_peak_sites)  # output the sites
            m_original_sites.clear()  #release the memory
            # YW added the following message
            print("Finished merging nearby clipped sites!")
            b_force = True
            rcd, basic_rcd = automatic_gnrt_parameters(sf_bam_list, sf_ref, s_working_folder, n_jobs,
                                                       b_force, b_tumor, f_purity)
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
            if b_automatic==True:
                n_disc_cutoff=rcd[1]
                # if b_tumor==True:
                #     n_disc_cutoff = adjust_cutoff_tumor(n_disc_cutoff, 0)
            print("Discordant cutoff: {0} is used!!!".format(n_disc_cutoff))

            sf_tmp = s_working_folder + "disc_tmp.list"
            # YW 2020/08/03 github update: added the following line
            sf_raw_disc=sf_out + global_values.RAW_DISC_TMP_SUFFIX #save the left and right raw disc for each site
            tem_locator = TE_Multi_Locator(sf_bam_list, s_working_folder, n_jobs, sf_ref)
            # YW 2020/08/03 github update: added 2 arguments sf_raw_disc, b_tumor
            tem_locator.filter_candidate_sites_by_discordant_pairs_multi_alignmts(m_sites_clip_peak, iextend, i_is, f_dev,
                                                                                  n_disc_cutoff, sf_annotation_Alu, sf_annotation_L1,
                                                                                  sf_annotation_SVA, sf_tmp, sf_raw_disc, b_tumor)
            # YW 2020/04/23 added the if statement to alert in advance the problem of unable to merge_clip_disc
            if os.stat(sf_tmp).st_size == 0:
                print("{0} is EMPTY! It will be impossible to merge clip disc in the next step. Set all counts of left and right discordant reads to 0. Discordant cutoff {1} is too high.\n".format(sf_tmp, n_disc_cutoff))
            # YW 2020/07/20 modified merge_clip_disc function to make sure locations without discordant read support will go through
            xfilter.merge_clip_disc(sf_tmp, sf_candidate_list, sf_out)
####
####
    ####
    elif args.filter_csn:  #filter out the FP by the pattern in the consensus repeat
        print("Working on \"clip-disc-filtering\" step!")
        sf_bam_list = args.bam  ###read in a bam list file
        s_working_folder = args.wfolder
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
            rcd, basic_rcd = automatic_gnrt_parameters(sf_bam_list, sf_ref, s_working_folder, n_jobs, b_force,
                                                       b_tumor, f_purity)
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
                # if b_tumor==True:
                #     n_clip_cutoff = adjust_cutoff_tumor(n_clip_cutoff)
                #     n_disc_cutoff = adjust_cutoff_tumor(n_disc_cutoff, 2)
            print("Filter (on cns) cutoff: number of clipped: {0} and number of discordant reads: {1} are used!!!\n".format(n_clip_cutoff, n_disc_cutoff))

            x_cd_filter = XClipDiscFilter(sf_bam_list, s_working_folder, n_jobs, sf_ref)
            x_cd_filter.call_MEIs_consensus(sf_candidate_list, sf_raw_disc, iextnd, bin_size,
                                            sf_rep_cns_Alu, sf_rep_cns_L1, sf_rep_cns_SVA,
                                            sf_flank, i_flank_lenth,
                                            bmapped_cutoff, i_concord_dist, f_concord_ratio, n_clip_cutoff, n_disc_cutoff,
                                            sf_output)
    ####
####
    #why a new sub-module?
    #1. we need to collect clip and disc reads for all the candidates (filter out some existing and not qualified ones)
    #2. The old module doesn't work well!
    #3. we need to select the candidate
    elif args.transduction:#
        #need to re-collect all the clip, disc reads
        sf_bam_list = args.bam  ###read in a bam list file
        sf_candidate_list = args.input #this is the output from the "cns" step.
        sf_raw_disc = args.input2  # this is the raw disc file
        iextnd = 400  ###for each site, re-collect reads in range [-iextnd, iextnd], this around ins +- 3*derivation
        bin_size = 50000000  # block size for parallelization
        sf_cns = args.reference  ####repeat copies/cns here
        s_working_folder = args.wfolder
        print("Current working folder is: {0}\n".format(s_working_folder))
        n_jobs = args.cores
        sf_reference = args.ref  ###reference genome, some cram file require this file to open
        sf_flank = args.fflank  # this is the flanking region
        i_flank_lenth = args.flklen
        sf_output = args.output
        sf_rmsk = args.annotation
        i_rep_type = args.rep_type
####
        if b_resume == False or os.path.isfile(sf_output) == False:
            if os.path.isfile(sf_flank)==True:#for Alu and many others, there is no transduction
                b_force = False
                rcd, basic_rcd = automatic_gnrt_parameters(sf_bam_list, sf_reference, s_working_folder, n_jobs, b_force,
                                                           b_tumor, f_purity)
                ave_cov = basic_rcd[0]  # ave coverage
                rlth = basic_rcd[1]  # read length
                mean_is = basic_rcd[2]  # mean insert size
                std_var = basic_rcd[3]  # standard derivation
                print("Mean insert size is: {0}\n".format(mean_is))
                print("Standard derivation is: {0}\n".format(std_var))
                max_is = int(mean_is + 3 * std_var)
                if iextnd < max_is:  # correct the bias
                    iextnd = max_is
                i_concord_dist=550
                f_concord_ratio = 0.25
                if i_concord_dist < max_is:  # correct the bias
                    i_concord_dist = max_is
                global_values.set_read_length(rlth)
                global_values.set_insert_size(max_is)
                global_values.set_average_cov(ave_cov)
                n_clip_cutoff = args.cliprep  # this is the sum of left and right clipped reads
                n_disc_cutoff = args.ndisc  # each sample should have at least this number of discordant reads
                if b_automatic == True:
                    n_clip_cutoff = rcd[0]
                    n_disc_cutoff = rcd[1]
                xtransduction = XTransduction(s_working_folder, n_jobs, sf_reference)
                i_win=150 #if a site is close to an existing site, then will not be considered again
                sf_tmp_slct=sf_raw_disc+".slct"
                i_min_copy_len = 225
                xannotation = xtransduction.prep_annotation_interval_tree(sf_rmsk, i_min_copy_len)
                xtransduction.re_slct_with_clip_raw_disc_sites(sf_raw_disc, sf_candidate_list, n_disc_cutoff, xannotation,
                                             i_rep_type, b_tumor, sf_tmp_slct)
                #now for the selected sites, re-evaluate each one
                x_cd_filter = XClipDiscFilter(sf_bam_list, s_working_folder, n_jobs, sf_reference)
                i_max_cov=ave_cov*(global_values.MAX_COV_TIMES+1)
                sf_output_tmp=sf_output + global_values.TD_NON_SIBLING_SUFFIX
                xtransduction.call_candidate_transduction_v3(sf_tmp_slct, sf_candidate_list, x_cd_filter,
                                                             sf_flank, sf_cns, i_flank_lenth, iextnd, bin_size, n_clip_cutoff,
                                                             n_disc_cutoff, i_concord_dist, f_concord_ratio, xannotation,
                                                             sf_bam_list, i_rep_type, i_max_cov, ave_cov, sf_output_tmp)
####
                xorphan=XOrphanTransduction(s_working_folder, n_jobs, sf_reference)
                n_half_disc_cutoff=n_disc_cutoff/2
                i_search_win=2000
                sf_updated_cns=sf_output #this is the final updated

                #1.Call out the sibling transduction events from the current list
                sf_sibling_TD=sf_output+".sibling_transduction_from_existing_list"
                xorphan.call_sibling_TD_from_existing_list(sf_output_tmp, sf_bam_list, iextnd, n_half_disc_cutoff,
                                                           i_search_win, xannotation, i_rep_type, i_max_cov,
                                                           sf_updated_cns, sf_sibling_TD)

                # #2. Call orphan "sibling" transdcution from non_existing list
                # sf_sibling_TD2 = sf_output + ".novel_sibling_transduction"
                # b_with_original=False
                # sf_tmp_slct2=sf_raw_disc+".slct2"
                # #select the sites to exclude the already called out sites
                # xorpha.re_slct_with_clip_raw_disc_sites(sf_raw_disc, sf_output_tmp, n_disc_cutoff, xannotation,
                #                                                i_rep_type, b_tumor, sf_tmp_slct2, b_with_original)
####
                #update high confident ones (in "cns" filter step, two results are generated)
                sf_ori_hc=sf_candidate_list+global_values.HIGH_CONFIDENT_SUFFIX
                sf_new_hc=sf_output+global_values.HIGH_CONFIDENT_SUFFIX
                xorphan.update_high_confident_callset(sf_ori_hc, sf_updated_cns, sf_new_hc)
####
            else:#rename the two files generated in previous step
                copyfile(sf_candidate_list, sf_output)
                sf_ori_hc = sf_candidate_list + global_values.HIGH_CONFIDENT_SUFFIX
                sf_new_hc=sf_output+global_values.HIGH_CONFIDENT_SUFFIX
                copyfile(sf_ori_hc, sf_new_hc)

####
    elif args.sibling:#sibling orphan transduction
        '''
        Todo: 09-29-2019: Add filtering modules:
        1. using background low mapq reads (multiple mapped reads) for filtering
        2. Set a upper-bound cutoff for discordant reads
        3. using blacklist for filtering
        '''
        sf_bam_list = args.bam  ###read in a bam list file
        sf_pre_step_out = args.input  # this is the output from the "cns" step.
        sf_raw_disc = args.input2  # this is the raw disc file
        iextnd = 400  ###for each site, re-collect reads in range [-iextnd, iextnd], this around ins +- 3*derivation
        bin_size = 50000000  # block size for parallelization
        sf_cns = args.reference  ####repeat copies/cns here
        s_working_folder = args.wfolder
        print("Current working folder is: {0}\n".format(s_working_folder))
        n_jobs = args.cores
        sf_reference = args.ref  ###reference genome, some cram file require this file to open
        sf_flank = args.fflank  # this is the flanking region
        i_flank_lenth = args.flklen
        sf_output = args.output
        sf_rmsk = args.annotation
        i_rep_type = args.rep_type
        sf_black_list = args.blacklist

        if os.path.isfile(sf_flank) == True:  #for Alu and many others, there is no transduction
            b_force = False
            rcd, basic_rcd = automatic_gnrt_parameters(sf_bam_list, sf_reference, s_working_folder, n_jobs, b_force,
                                                       b_tumor, f_purity)
            ave_cov = basic_rcd[0]  # ave coverage
            rlth = basic_rcd[1]  # read length
            mean_is = basic_rcd[2]  # mean insert size
            std_var = basic_rcd[3]  # standard derivation
            print("Mean insert size is: {0}\n".format(mean_is))
            print("Standard derivation is: {0}\n".format(std_var))
            max_is = int(mean_is + 3 * std_var)

            i_concord_dist = 550
            f_concord_ratio = 0.25
            if i_concord_dist < max_is:  # correct the bias
                i_concord_dist = max_is
            global_values.set_read_length(rlth)
            global_values.set_insert_size(mean_is) #here set mean inset size
            global_values.set_average_cov(ave_cov)

            n_clip_cutoff = args.cliprep  # this is the sum of left and right clipped reads
            n_disc_cutoff = args.ndisc  # each sample should have at least this number of discordant reads
            if b_automatic == True:
                n_clip_cutoff = rcd[0]
                n_disc_cutoff = rcd[1]

            xorphan = XOrphanTransduction(s_working_folder, n_jobs, sf_reference)
            i_min_copy_len = 225
            xorphan.set_boundary_extend(mean_is)
            xannotation = xorphan.prep_annotation_interval_tree(sf_rmsk, i_min_copy_len)
            # 2. Call orphan "sibling" transdcution from non_existing list
            # sf_sibling_TD2 = sf_output + ".novel_sibling_transduction"
            b_with_original = False
            sf_tmp_slct2 = sf_raw_disc + ".slct2"
            sf_output_tmp = sf_pre_step_out + global_values.TD_NON_SIBLING_SUFFIX
            # select the sites to exclude the already called out sites, and filter out sites fall in black_list
            xorphan.re_slct_with_clip_raw_disc_sites(sf_raw_disc, sf_output_tmp, n_disc_cutoff, xannotation,
                                                    i_rep_type, b_tumor, sf_tmp_slct2, b_with_original)

            #re-select transduction candidates based on disc-clip consistency (clip position encompass disc ones?)
            m_failed_ori_td=xorphan.distinguish_source_from_insertion_for_td(sf_pre_step_out, sf_bam_list, i_concord_dist,
                                                                           n_clip_cutoff, n_disc_cutoff, sf_black_list,
                                                                           sf_rmsk)
            sf_sibling_TD=sf_output
            if os.path.isfile(sf_black_list)==False:
                print("Blacklist file {0} does not exist!".format(sf_black_list))
            xorphan.call_novel_sibling_TD_from_raw_list(sf_tmp_slct2, sf_bam_list, i_concord_dist, n_clip_cutoff,
                                                        n_disc_cutoff, sf_black_list, sf_rmsk, sf_sibling_TD)

            ####append the newly called events to existing list
            xorphan.append_to_existing_list(sf_sibling_TD, sf_pre_step_out, m_failed_ori_td)
            xorphan.append_to_existing_list(sf_sibling_TD, sf_pre_step_out+global_values.HIGH_CONFIDENT_SUFFIX,
                                            m_failed_ori_td)
        else:
            #sf_output_tmp = sf_pre_step_out + global_values.TD_NON_SIBLING_SUFFIX
            if os.path.isfile(sf_pre_step_out)==True:#do td filtering only
                b_force = False
                rcd, basic_rcd = automatic_gnrt_parameters(sf_bam_list, sf_reference, s_working_folder, n_jobs, b_force,
                                                           b_tumor, f_purity)
                ave_cov = basic_rcd[0]  # ave coverage
                rlth = basic_rcd[1]  # read length
                mean_is = basic_rcd[2]  # mean insert size
                std_var = basic_rcd[3]  # standard derivation
                print("Mean insert size is: {0}\n".format(mean_is))
                print("Standard derivation is: {0}\n".format(std_var))
                max_is = int(mean_is + 3 * std_var)

                i_concord_dist = 550
                f_concord_ratio = 0.25
                if i_concord_dist < max_is:  # correct the bias
                    i_concord_dist = max_is
                global_values.set_read_length(rlth)
                global_values.set_insert_size(mean_is)  # here set mean inset size
                global_values.set_average_cov(ave_cov)
                n_clip_cutoff = args.cliprep  # this is the sum of left and right clipped reads
                n_disc_cutoff = args.ndisc  # each sample should have at least this number of discordant reads
                if b_automatic == True:
                    n_clip_cutoff = rcd[0]
                    n_disc_cutoff = rcd[1]
                xorphan = XOrphanTransduction(s_working_folder, n_jobs, sf_reference)
                # re-select transduction candidates based on disc-clip consistency (clip position encompass disc ones?)
                m_failed_ori_td = xorphan.distinguish_source_from_insertion_for_td(sf_pre_step_out, sf_bam_list,
                                                                                   i_concord_dist, n_clip_cutoff,
                                                                                   n_disc_cutoff, sf_black_list, sf_rmsk)
                xorphan.update_existing_list_only(sf_pre_step_out, m_failed_ori_td)
####
####
    ####this module for: 1) tumor case-control files;
    ####2) trio or quads to call de novo insertion
    elif args.case_control:#case-control mode to call somatic events
        b_somatic_hc = args.somatic_hc
        sf_candidate_list = args.input  # this is the list called from case
        sf_output = args.output
        s_working_folder = args.wfolder
        if b_somatic_hc==False:#this is for all the raw call set
            sf_bam_list = args.bam  # this is the control bam file list
            sf_ref = args.ref  # reference genome
            n_jobs = args.cores
            nclip_cutoff = args.cliprep  # this is the sum of left and right clipped reads
            ndisc_cutoff = args.ndisc  # each sample should have at least this number of discordant reads
            sf_rep_cns = args.reference  ####repeat copies/cns here
            sf_flank = args.fflank  # this is the flanking region
            i_flk_len = args.flklen
            b_force=True
            rcd=None
            basic_rcd=None
            n_polyA_cutoff=0
            if b_automatic==True:
                rcd, basic_rcd=automatic_gnrt_parameters_case_control(sf_bam_list, sf_ref, s_working_folder, n_jobs, b_force)
                nclip_cutoff=rcd[0]
                ndisc_cutoff=rcd[1]
                n_polyA_cutoff=rcd[2]
            ccm=CaseControlMode(sf_ref, s_working_folder, n_jobs)
            #ccm.set_parameters(iextnd, bin_size, bmapped_cutoff, i_concord_dist, f_concord_ratio)
            rlth = basic_rcd[1]  # read length
            mean_is = basic_rcd[2]  # mean insert size
            std_var = basic_rcd[3]  # standard derivation
            max_is = int(mean_is + 3 * std_var) + int(rlth)
            extnd = max_is
            bin_size = 50000000  # block size for parallelization
            print("clip,disc,polyA-cutoff is ({0}, {1}, {2})".format(nclip_cutoff, ndisc_cutoff, n_polyA_cutoff))
            n_polyA_cutoff=ndisc_cutoff #if both sides have more than cutoff polyA, then filter out
            ccm.call_somatic_TE_insertion(sf_bam_list, sf_candidate_list, extnd, nclip_cutoff, ndisc_cutoff,
                                          n_polyA_cutoff, sf_rep_cns, sf_flank, i_flk_len, bin_size, sf_output, b_tumor)
    ####
        else:#This is to parse out the high confident somatic ones, assume already have the raw somatic callset
            sf_raw_somatic=args.input2
            sf_ref=""
            n_jobs=1
            ccm = CaseControlMode(sf_ref, s_working_folder, n_jobs)
            ccm.parse_high_confident_somatic(sf_candidate_list, sf_raw_somatic, sf_output)

####
    elif args.mosaic:  # this is only for normal illumina data
        #for mosaic events, when check clip information, we will check the polyA information
        print("Working on mosaic \"clip\" step!")
        sf_bam_list = args.input
        s_working_folder = args.wfolder
        n_jobs = args.cores
        sf_rep_cns = args.cns
        sf_rep = args.reference  ####repeat copies "-r"
        sf_annotation = args.annotation
        sf_out = args.output
        b_se = args.single  ##single end reads or not, default is not
        sf_ref = args.ref  ###reference genome "-ref"
        b_force = args.force  # force to run from the very beginning
        if b_force == True:
            global_values.set_force_clean()
        site_clip_cutoff = args.siteclip  # this is the cutoff for the exact position, use larger value for 10X
        global_values.set_initial_min_clip_cutoff(site_clip_cutoff)

        # merge the list from different bams of the same individual
        # Here when do the filtering, nearby regions are already considered!
        cutoff_left_clip = args.lclip
        cutoff_right_clip = args.rclip
        cutoff_clip_mate_in_rep = args.cliprep
        cutoff_polyA=1

        if b_automatic == True:
            rcd, basic_rcd = automatic_gnrt_parameters(sf_bam_list, sf_ref, s_working_folder, n_jobs)
            cutoff_left_clip = rcd[0]
            cutoff_right_clip = rcd[0]
            cutoff_clip_mate_in_rep = rcd[2]
        print("Clip cutoff: {0}, {1}, {2} are used!!!".format(cutoff_left_clip, cutoff_right_clip,
                                                              cutoff_clip_mate_in_rep))
        tem_locator = TE_Multi_Locator(sf_bam_list, s_working_folder, n_jobs, sf_ref)

        ####by default, if number of clipped reads is larger than this value, then discard
        max_cov_cutoff = int(0.5 * args.cov)  # at most half of coverage of clipped reads
        wfolder_pub_clip = args.cwfolder  # public clip folder
        ##Hard code inside:
        # 1. call_TEI_candidate_sites_from_clip_reads_v2 --> run_cnt_clip_part_aligned_to_rep_by_chrm_sort_version
        # here if half of the seq is mapped, then consider it as aligned work.
        ##2. require >=2 clip reads, whose clipped part is aligned to repeat copies
        tem_locator.call_TEI_candidate_sites_from_multiple_alignmts_mosaic(sf_annotation, sf_rep_cns, sf_rep, b_se,
                                                                    cutoff_left_clip,
                                                                    cutoff_right_clip, cutoff_clip_mate_in_rep,
                                                                    cutoff_polyA, wfolder_pub_clip,
                                                                    b_force, max_cov_cutoff, sf_out)
####
    elif args.barcode:  # this is only for 10X alignmt
        sf_ori_bam = args.bam  # pos indexted bam
        b_barcode = args.barcode
        sf_barcode_bam = args.barcode_bam  # barcode indexed bam
        s_working_folder = args.wfolder
        n_jobs = args.cores
        sf_annotation = args.annotation
        sf_candidate_list = args.input
        sf_out = args.output
        sf_ref = args.ref  ###reference genome, some cram file require this file to open

        xfilter = XIntemediateSites()
        m_input_sites = xfilter.load_in_candidate_list(sf_candidate_list)

        iextend = 2500
        i_cov_cutoff=automatic_set_molecule_cutoff_for_10X_bam(sf_ori_bam, sf_ref, s_working_folder, n_jobs)
        print("Using molecule coverage cutoff: {0}".format(i_cov_cutoff)) 
        #i_cov_cutoff = args.nbarcode ####here use an automatic one

        if s_working_folder[-1] != "/":
            s_working_folder += "/"
        sf_tmp = s_working_folder + "barcode_tmp.list"
        caller = TELocator(sf_ori_bam, sf_barcode_bam, s_working_folder, n_jobs, sf_ref)
        m_sites_barcode = caller.filter_candidate_sites_by_barcode_coverage(m_input_sites, iextend, i_cov_cutoff)
        caller.output_candidate_sites(m_sites_barcode, sf_tmp)  # clip_peak_discord_candidate_barcode.list
        sf_tmp2 = s_working_folder + "barcode_tmp2.list"
        xfilter.merge_clip_disc_barcode(sf_tmp, sf_candidate_list, sf_tmp2)

        ###combine the ones are close to each other
        window_size = 200  # if closer than 200, then the sites are merged as one
        xfilter.combine_closing_sites(sf_tmp2, window_size, sf_out)

####
    elif args.postF:#post filtering step
        s_working_folder = args.wfolder
        n_jobs = args.cores
        sf_xtea_rslt = args.input
        sf_rmsk = args.annotation
        #sf_rmsk_all=args.annotation2
        sf_new_out = args.output
        i_rep_type=args.rep_type
        sf_black_list=args.blacklist

        i_min_copy_len=225 #when check whether fall in repeat region, require the minimum copy length
        b_pf_mosaic=args.postFmosaic
        if b_pf_mosaic is True:#for mosaic events
            xpf_mosic = MosaicCaller(s_working_folder, n_jobs)
            xpf_mosic.run_call_mosaic(sf_xtea_rslt, sf_rmsk, i_min_copy_len, i_rep_type, sf_black_list, sf_new_out)
        else:
            xbasic=X_BasicInfo(s_working_folder)
            f_cov, f_ave_cov=xbasic.load_get_cov()
            xpost_filter = XPostFilter(s_working_folder, n_jobs)
            #here sf_black_list is the centromere + duplication region
            xpost_filter.run_post_filtering(sf_xtea_rslt, sf_rmsk, i_min_copy_len, i_rep_type, f_cov, sf_black_list,
                                            sf_new_out, b_tumor)
####
    ####
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
            sf_00_list = "/n/data1/hms/dbmi/park/simon_chu/projects/XTEA/genotyping/training_set_SSC/Genotyping/rslt_list/all_00.list"
            sf_01_list = "/n/data1/hms/dbmi/park/simon_chu/projects/XTEA/genotyping/training_set_SSC/Genotyping/rslt_list/all_01.list"
            sf_11_list = "/n/data1/hms/dbmi/park/simon_chu/projects/XTEA/genotyping/training_set_SSC/Genotyping/rslt_list/all_11.list"
            sf_arff = "/n/data1/hms/dbmi/park/simon_chu/projects/XTEA/genotyping/training_set_SSC/Genotyping/merged_all_0_1_2.arff"
            gc = GntpClassifier()
            #pkl_filename = "./genotyping/trained_model_ssc_py2_random_forest_two_category.pkl"
            gc.train_model(sf_arff, sf_model, f_ratio=0.01)
        else:#predict the genotype
            #sf_model = "./genotyping/trained_model_ssc_py2_random_forest_two_category.pkl"
            sf_xTEA = args.input #input raw results before calling genotype
            sf_new = args.output
            gc = GntpClassifier()
            pkl_model = gc.load_model_from_file(sf_model)
            sf_arff = sf_xTEA + ".arff"
            gc.predict_for_site(pkl_model, sf_xTEA, sf_new)

####
    elif args.gVCF:
        sf_ref = args.ref
        sf_bam_list = args.bam  ###read in a bam list file
        sf_raw_rslt = args.input
        i_rep_type=args.rep_type##
        sf_prefix=args.output
        sf_bam=""
        s_sample_id="null"

        with open(sf_bam_list) as fin_blist:
            for line in fin_blist:
                fields=line.split()
                if len(fields)<=0:
                    break
                sf_bam=fields[0]
                bam_fields=sf_bam.split("/")
                if len(bam_fields)<=0:
                    break
                s_sample_id2=bam_fields[-1]
                tmp_fields=s_sample_id2.split(".")
                s_sample_id=".".join(tmp_fields[:-1])
                break
        if os.path.isfile(sf_bam)==True:
            gvcf=gVCF()
            if sf_prefix[-1]!="/":
                sf_prefix+="/"
            reptype = RepType()
            s_rep_type = reptype.get_rep_type(i_rep_type)
            sf_vcf=sf_prefix+s_sample_id+"_"+s_rep_type+".vcf"
            gvcf.cvt_raw_rslt_to_gvcf(s_sample_id, sf_bam, sf_raw_rslt, i_rep_type, sf_ref, sf_vcf)
        else:
            print("Wrong bam file: {0}".format(sf_bam))
####
    elif args.joint:
        s_working_folder = args.wfolder
        n_jobs = args.cores
        sf_rslt_list = args.input
        sf_out=args.output
        islack=50
        if args.postFmosaic==True:
            mzc_joint=MosaicJointCalling(s_working_folder, n_jobs)
            mzc_joint.call_mosaic_from_multi_samples(sf_rslt_list, islack, sf_out)
####
    elif args.igv:#prepare the igv screenshot script for multiple bams and sites
        sf_sites = args.input #site list
        sf_bam_list=args.bam #bam list
        s_screenshot_folder = args.wfolder
        sf_out = args.output #output file
        sf_gnm=args.ref #"hg19" or "hg38"
        i_extnd = args.extend #"-e", "--extend"
        b_single_sample=args.single_sample
        x_igv = XIGV()
        if b_single_sample==True:
            sf_sites=x_igv.gnrt_sites_single_sample(sf_sites, sf_bam_list)
        x_igv.prepare_igv_scripts_multi_bams(sf_sites, sf_bam_list, s_screenshot_folder, i_extnd, sf_gnm, sf_out)

####
    elif args.collect:  # collect the reads for each candidate site
        sf_ori_bam = args.bam
        sf_barcode_bam = args.barcode_bam
        s_working_folder = args.wfolder
        n_jobs = args.cores
        sf_candidate_list = args.input
        i_cov_cutoff = args.nbarcode
        sf_annotation = args.annotation
        i_flank_lenth = args.flklen #extend the boundary of each annotated repeat
        sf_ref = args.ref  ###reference genome, some cram file require this file to open

        # collect reads for all the sites
        i_extend = 1500 #by default, collect barcode in [-1500, 1500] region
        xread_collection=XReadsCollection(s_working_folder, sf_ref)
        xread_collection.collect_phased_reads_all_TEIs(sf_ori_bam, sf_barcode_bam, sf_candidate_list, i_extend,
                                                       i_cov_cutoff, sf_annotation, i_flank_lenth, n_jobs)
####
    elif args.mutation:#call out the internal mutations by aligning the reads
        s_working_folder = args.wfolder
        n_jobs = args.cores
        sf_cns = args.reference  ####repeat copies/cns here
        sf_sites=args.input
        sf_merged_vcf=args.output
        n_len_cutoff=args.teilen
        xmutation=XMutation(s_working_folder)
        xmutation.call_mutations_from_reads_algnmt(sf_sites, sf_cns, n_len_cutoff, n_jobs, sf_merged_vcf)

####
    elif args.gene:#
        #To run: --gene -a -i -o
        sf_gene_annotation=args.annotation
        sf_input=args.input
        sf_output=args.output
        gff=GFF3(sf_gene_annotation)
        iextnd=global_values.UP_DOWN_GENE
        gff.load_gene_annotation_with_extnd(iextnd)
        gff.index_gene_annotation_interval_tree()
        gff.annotate_results(sf_input, sf_output)

    # this step needs to retrieve the seqs of mate reads
    elif args.collect_illumina:  # this is only for normal illumina data
        sf_ori_bam = args.bam
        sf_candidate_list = args.input
        n_jobs = args.cores
        s_working_folder = args.wfolder

    elif args.assembly:  # assemble the reads for all the sites
        b_local = args.local
        s_working_folder = args.wfolder
        n_jobs = args.cores
        sf_candidate_list = args.input
        sf_ref = args.ref  ###reference genome, some cram file require this file to open

        xlasm = XLocalAssembly(s_working_folder, sf_ref)
        if b_local == True:
            # xlasm.assemble_all_TEIs_locally(sf_candidate_list, n_jobs)
            xlasm.assemble_all_phased_TEIs_locally(sf_candidate_list, n_jobs)
        else:
            # xlasm.assemble_all_TEIs_slurm_cluster(sf_candidate_list)
            xlasm.assemble_all_phased_TEIs_slurm_cluster(sf_candidate_list)
####
    elif args.map:####align the asm to reference genome
        sf_ref = args.reference#
        s_working_folder = args.wfolder
        n_jobs = int(args.cores)
        sf_sites = args.input
        sf_cns=args.ref####repeat consensus
        sf_final_list=args.output
        sf_final_seqs=sf_final_list+".fa"

        xctg = XTEContig(s_working_folder, n_jobs)
        xctg.align_asm_contigs_to_reference(sf_sites, sf_ref, s_working_folder)
        xctg.call_MEIs_from_all_group_contigs(sf_sites, sf_ref, sf_cns, s_working_folder, sf_final_list, sf_final_seqs)

    elif args.flk_map:  #gnrt the flank regions and align the regions to the contigs 
        sf_ref = args.ref
        s_working_folder = args.wfolder
        n_jobs = int(args.cores)
        sf_sites = args.input

        i_extend = 500
        xref = XReference()
        xref.gnrt_flank_region_for_sites(sf_sites, i_extend, n_jobs, sf_ref, s_working_folder)
        xctg = XTEContig(s_working_folder, n_jobs)
        xctg.align_flanks_to_phased_contig(sf_sites)
####
    elif args.filter_asm:
        s_working_folder = args.wfolder
        n_jobs = int(args.cores)
        sf_sites = args.input
        sf_keep_sites = args.output
        sf_repeat_copies = args.reference  ####repeat copies here

        # xctg.filter_out_non_TE_from_asm(sf_sites, flank_length, f_map_cutoff, i_slack, sf_keep_sites)
        flank_length = 500
        f_map_cutoff = 0.45 #25% of the bases are required to be mapped! This is a small ratio
        i_slack = 35
        xctg = XTEContig(s_working_folder, n_jobs)
        xctg.validate_TEI_from_phased_asm_algnmt(sf_sites, flank_length, f_map_cutoff, i_slack, sf_repeat_copies,
                                                 sf_keep_sites)
    elif args.sv:
        sf_sites = args.input
        sf_dels=args.output
        te_del=TEDeletion()
        m_del=te_del.call_del_matched_brkpnts(sf_sites)
        te_del.output_del(m_del, sf_dels)

    elif args.collect_clip:#collect the clipped reads for the sample
        sf_bam_list = args.input
        s_working_folder = args.wfolder ##this is the folder to save all the clipped reads of the sample
        n_jobs = args.cores
        b_se = args.single  ##single end reads or not, default is not
        sf_ref = args.ref  ###reference genome "-ref"
        sf_annotation = args.annotation

        tem_locator = TE_Multi_Locator(sf_bam_list, s_working_folder, n_jobs, sf_ref)
        s_clip_wfolder=s_working_folder
        wfolder_pub_clip = args.cwfolder  # public clip folder
        # collect the clipped reads only
        tem_locator.collect_all_clipped_from_multiple_alignmts(sf_annotation, b_se, s_clip_wfolder, wfolder_pub_clip)

    elif args.visualization:  ##show the shared barcode heatmap between the TEI site and source region
        sf_matrix = args.input
        ####classify the TE insertions to normal ones and complex ones:
        # for insertion with deletion: the other clipped part will not be aligned to a repeat region, but the other breakpoint
        # while for normal ones, the other clipped part will be aligned to a repeat region
        ####2. python x_TEA_main.py -T -i ${VCF} -r ${ANNOTATION} -b ${BAM} -d ${BARCODE_BAM} -p ${TMP} -n 1
        ####4.
####
