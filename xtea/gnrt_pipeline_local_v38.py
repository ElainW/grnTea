#!/usr/bin/env python3

##11/04/2018
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong.simon.chu@gmail.com

import os
import sys
from subprocess import *
import argparse
import ntpath
import global_values
ILLUMINA="illumina"
X10="10X"
S_VERSION="v0.1"

####
PUB_CLIP="pub_clip"
PUB_CLIP_LOCUS="pub_clip_locus"
REP_TYPE_L1="L1"
REP_TYPE_ALU="Alu"
REP_TYPE_SVA="SVA"
REP_TYPE_HERV="HERV"
REP_TYPE_MIT="Mitochondrion"
REP_TYPE_MSTA="MSTA"
CASE_LIST_SUFFIX=".case_list"
CTRL_LIST_SUFFIX=".ctrl_list"

# YW 2020/07/24 customized memory usage by rep_type
# YW 2021/05/20 removed rep_type
def gnrt_script_head(spartition, ncores, stime, smemory, s_id):
    s_head = "#!/bin/bash\n\n"
    s_head += f"#SBATCH -c {ncores}\n"
    s_head += f"#SBATCH -t {stime}\n"
    # if rep_type == REP_TYPE_L1:
    #     memory = 10
    #     # memory = 20 # for highcov
    # elif rep_type == REP_TYPE_SVA:
    #     memory = 5
    #     # memory = 15 # for highcov
    # else:
    #     memory = smemory
    s_head += f"#SBATCH --mem={smemory}G\n"
    s_head += f"#SBATCH -p {spartition}\n"
    s_head += f"#SBATCH -J {s_id}\n" # YW added 2021/05/24
    s_head += f"#SBATCH -o {s_id}/%j.out\n"
    s_head += "#SBATCH --mail-type=END,FAIL\n"
    s_head += "#SBATCH --mail-user=yilanwang@g.harvard.edu\n"
    if spartition == "park" or spartition == "priopark":
        s_head += "#SBATCH --account=park_contrib\n\n"
    return s_head

####
def gnrt_script_head_lsf(spartition, ncores, stime, smemory, s_id):
    s_head = "#!/bin/bash\n\n"
    s_head += f"#BSUB -n {ncores}\n"
    if stime!=None and stime!="":
        s_head += f"#BSUB -W {stime}\n"
    # if rep_type == REP_TYPE_L1:
    #     memory = 10
    # elif rep_type == REP_TYPE_SVA:
    #     memory = 5
    # else:
    #     memory = smemory
    s_head += f"#BSUB -M {memory}G\n"
    s_head += f"#BSUB -q {spartition}\n"
    s_head += f"#BSUB -o {s_id}_%J.out\n"
    s_head += f"#BSUB -e {s_id}_%J.err\n"
    s_head += f"#BSUB -R \"rusage[mem={int(memory)/int(ncores)}G]\"\n" #this is memory per core
    return s_head


# load in the parameter file or the configuration file
def load_par_config(sf_par_config):
    # by default, SF_FLANK is set to null, as Alu no need for SF_FLANK, as we don't check transduction for Alu
    l_pars = []
    with open(sf_par_config) as fin_par_config:
        for line in fin_par_config:
            if len(line) > 0 and line[0] == "#":
                continue
            fields = line.split()
            l_pars.append((fields[0], fields[1]))
    return l_pars

# gnrt pars
def gnrt_parameters(l_pars):
    s_pars = ""
    for rcd in l_pars:
        sid = rcd[0]
        svalue = str(rcd[1])
        sline = sid + "=" + svalue + "\n"
        s_pars += sline
    return s_pars
####
####
# grnt calling steps
# YW 2021/05/20 added icns_c, took out b_tumor, f_purity, b_mosaic, changed i_rep_type to l_rep_type
# YW 2021/05/23 changed sclip_step, sdisc_step command, add annotation/reference/cns by rep_type in l_rep_type (Alu/L1/SVA), commented out other steps and excluded other rep type (raise NonImplementedError)
# YW 2021/08/04 added b_ctrl, sf_ref_bed, error_margin to enable control bam file coordinate lifting
def gnrt_calling_command(iclip_c, iclip_rp, idisc_c, icns_c, ncores, iflk_len, iflag,
                         b_user_par, b_force, b_resume, l_rep_type, s_cfolder, s_cfolder_locus,
                         c_realn_partition, c_realn_time, c_realn_mem, d_realn_partition, d_realn_time, d_realn_mem, check_interval,
                         b_ctrl, sf_ref_bed, error_margin, sf_locus_file, feat_extract_time, feat_extract_partition, email_user):
    s_user = ""
    if b_user_par == True:
        s_user = "--user"
    s_clean = ""
    if b_force == True:
        s_clean = "--force"
    # s_tumor = ""
    # s_purity=""
    # if b_tumor == True:
    #     s_tumor = "--tumor"
    #     s_purity="--purity {0}".format(f_purity)
    s_resume = ""
    if b_resume == True:
        s_resume = "--resume"
    
    if REP_TYPE_HERV in l_rep_type or REP_TYPE_MSTA in l_rep_type or REP_TYPE_MIT in l_rep_type:
        raise NotImplementedError
    setup_step = "cd ${PREFIX}\n"
    setup_step += "ln -nsf ${XTEA_PATH} ${TMP}source_scripts\n" # YW 2021/05/26 added to import modules from working directory
    sclip_step = f"python3 ${{XTEA_PATH}}\"x_TEA_main.py\" -C -i ${{BAM_LIST}} --lc {iclip_c} --rc {iclip_c} --cr {iclip_rp} " \
                 f"-o ${{PREFIX}}\"candidate_list_from_clip.txt\"  -n {ncores} --cp {s_cfolder} {s_user} {s_clean} {s_resume} " \
                 f"--ref ${{REF}} -p ${{TMP}} --c_realn_partition {c_realn_partition} --c_realn_time {c_realn_time} --c_realn_mem {c_realn_mem} --check_interval {check_interval} --email_user={email_user} "
    if REP_TYPE_ALU in l_rep_type:
        sclip_step += "--Alu-annotation ${ALU_ANNOTATION} --Alu-reference ${ALU_COPY_WITH_FLANK} --Alu-cns ${ALU_CNS} "
    if REP_TYPE_L1 in l_rep_type:
        sclip_step += "--L1-annotation ${L1_ANNOTATION} --L1-reference ${L1_COPY_WITH_FLANK} --L1-cns ${L1_CNS} "
    if REP_TYPE_SVA in l_rep_type:
        sclip_step += "--SVA-annotation ${SVA_ANNOTATION} --SVA-reference ${SVA_COPY_WITH_FLANK} --SVA-cns ${SVA_CNS}\n"
    # if b_SVA is True:
    #     sclip_step = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" -C --sva -i ${{BAM_LIST}} --lc {0} --rc {1} --cr {2}  " \
    #                  "-r ${{L1_COPY_WITH_FLANK}}  -a ${{ANNOTATION}} --cns ${{L1_CNS}} --ref ${{REF}} -p ${{TMP}} " \
    #                  "-o ${{PREFIX}}\"candidate_list_from_clip.txt\"  -n {3} --cp {4} {5} {6} {7} {8} {9}\n".format(iclip_c, iclip_c,
    #                                                                                                 iclip_rp, ncores,
    #                                                                                                 s_cfolder, s_purity, s_user, s_clean, s_tumor, s_resume)
    sdisc_step = f"python3 ${{XTEA_PATH}}\"x_TEA_main.py\" -D -i ${{PREFIX}}\"candidate_list_from_clip.txt\" --nd {idisc_c} " \
                 f"--ncns {icns_c} --ref ${{REF}} -b ${{BAM_LIST}} -p ${{TMP}} -o ${{PREFIX}}\"candidate_list_from_disc.txt\" -n {ncores} -m ${{PREFIX}}\"feature_matrix.txt\" {s_user} {s_resume} " \
                 f"--d_realn_partition {d_realn_partition} --d_realn_time {d_realn_time} --d_realn_mem {d_realn_mem} --check_interval {check_interval} --error_margin {error_margin} --email_user={email_user} " # YW 2021/08/04 added error_margin
    if REP_TYPE_ALU in l_rep_type:
        sdisc_step += "--Alu-annotation ${ALU_ANNOTATION} --Alu-cns ${ALU_CNS} "
    if REP_TYPE_L1 in l_rep_type:
        sdisc_step += "--L1-annotation ${L1_ANNOTATION} --L1-cns ${L1_CNS} "
    if REP_TYPE_SVA in l_rep_type:
        sdisc_step += "--SVA-annotation ${SVA_ANNOTATION} --SVA-cns ${SVA_CNS} "
    # YW 2021/08/04 added to enable control bam file coordinate lifting
    if args.train: # YW 2021/10/28 added args.train if statement
        if b_ctrl:
            sdisc_step += f"--train --ctrl --ref_bed {sf_ref_bed}\n"
        else:
            sdisc_step += "--train --ctrl_bed ${CTRL_BED}\n"
    else:
        sdisc_step += f"--feat_extract_time {feat_extract_time} --feat_extract_partition {feat_extract_partition}\n"
    
    sclip_locus_step = f"python3 ${{XTEA_PATH}}\"x_TEA_main.py\" -L --locus_file {sf_locus_file} -i ${{BAM_LIST}} --lc {iclip_c} --rc {iclip_c} --cr {iclip_rp} " \
                 f"-o ${{PREFIX}}\"candidate_list_from_clip_locus.txt\"  -n {ncores} --cp {s_cfolder_locus} {s_user} {s_clean} {s_resume} " \
                 f"--ref ${{REF}} -p ${{TMP}} --c_realn_partition {c_realn_partition} --c_realn_time {c_realn_time} --c_realn_mem {c_realn_mem} --check_interval {check_interval} --email_user={email_user} "
    if REP_TYPE_ALU in l_rep_type:
        sclip_locus_step += "--Alu-annotation ${ALU_ANNOTATION} --Alu-reference ${ALU_COPY_WITH_FLANK} --Alu-cns ${ALU_CNS} "
    if REP_TYPE_L1 in l_rep_type:
        sclip_locus_step += "--L1-annotation ${L1_ANNOTATION} --L1-reference ${L1_COPY_WITH_FLANK} --L1-cns ${L1_CNS} "
    if REP_TYPE_SVA in l_rep_type:
        sclip_locus_step += "--SVA-annotation ${SVA_ANNOTATION} --SVA-reference ${SVA_COPY_WITH_FLANK} --SVA-cns ${SVA_CNS}\n"
    
    sdisc_locus_step = f"python3 ${{XTEA_PATH}}\"x_TEA_main.py\" -D -i ${{PREFIX}}\"candidate_list_from_clip_locus.txt\" --nd {idisc_c} " \
                       f"--ncns {icns_c} --ref ${{REF}} -b ${{BAM_LIST}} -p ${{TMP}} -o ${{PREFIX}}\"candidate_list_from_disc_locus.txt\" -n {ncores} -m ${{PREFIX}}\"feature_matrix_locus.txt\" {s_user} {s_resume} " \
                       f"--d_realn_partition {d_realn_partition} --d_realn_time {d_realn_time} --d_realn_mem {d_realn_mem} --check_interval {check_interval} --error_margin {error_margin} --email_user={email_user} "
    if REP_TYPE_ALU in l_rep_type:
        sdisc_locus_step += "--Alu-annotation ${ALU_ANNOTATION} --Alu-cns ${ALU_CNS} "
    if REP_TYPE_L1 in l_rep_type:
        sdisc_locus_step += "--L1-annotation ${L1_ANNOTATION} --L1-cns ${L1_CNS} "
    if REP_TYPE_SVA in l_rep_type:
        sdisc_locus_step += "--SVA-annotation ${SVA_ANNOTATION} --SVA-cns ${SVA_CNS} "
    # YW 2021/08/04 added to enable control bam file coordinate lifting
    if args.train: # YW 2021/10/28 added args.train if statement
        if b_ctrl:
            sdisc_locus_step += f"--train --ctrl --ref_bed {sf_ref_bed}\n"
        else:
            sdisc_locus_step += "--train --ctrl_bed ${CTRL_BED}\n"
    # if b_SVA is True:
    #     sdisc_step = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\"  -D --sva -i ${{PREFIX}}\"candidate_list_from_clip.txt\" --nd {0} " \
    #                  "--ref ${{REF}} -a ${{ANNOTATION}} -b ${{BAM_LIST}} -p ${{TMP}} " \
    #                  "-o ${{PREFIX}}\"candidate_list_from_disc.txt\" -n {1} {2} {3} {4} {5}\n".format(idisc_c, ncores, s_purity, s_user, s_tumor, s_resume)
    # if b_mosaic==True:
    #     sclip_step = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" --mosaic -i ${{BAM_LIST}} --lc {0} --rc {1} --cr {2}  " \
    #                  "-r ${{L1_COPY_WITH_FLANK}}  -a ${{ANNOTATION}} --cns ${{L1_CNS}} --ref ${{REF}} -p ${{TMP}} " \
    #                  "-o ${{PREFIX}}\"candidate_list_from_clip.txt\"  -n {3} --cp {4} {5} {6}\n".format(iclip_c, iclip_c,
    #                                                                                                 iclip_rp,
    #                                                                                                 ncores, s_cfolder,
    #                                                                                                 s_user, s_clean)
    #     sdisc_step = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\"  -D --somatic -i ${{PREFIX}}\"candidate_list_from_clip.txt\" --nd {0} " \
    #                  "--ref ${{REF}} -a ${{ANNOTATION}} -b ${{BAM_LIST}} -p ${{TMP}} " \
    #                  "-o ${{PREFIX}}\"candidate_list_from_disc.txt\" -n {1} {2}\n".format(idisc_c, ncores, s_user)
    # if b_mosaic is True and b_SVA is True:
    #     sclip_step = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" --mosaic --sva -i ${{BAM_LIST}} --lc {0} --rc {1} --cr {2}  " \
    #                  "-r ${{L1_COPY_WITH_FLANK}}  -a ${{ANNOTATION}} --cns ${{L1_CNS}} --ref ${{REF}} -p ${{TMP}} " \
    #                  "-o ${{PREFIX}}\"candidate_list_from_clip.txt\"  -n {3} --cp {4} {5} {6}\n".format(iclip_c, iclip_c,
    #                                                                                                 iclip_rp, ncores,
    #                                                                                                 s_cfolder, s_user, s_clean)
    #     sdisc_step = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\"  -D --sva --somatic -i ${{PREFIX}}\"candidate_list_from_clip.txt\" --nd {0} " \
    #                  "--ref ${{REF}} -a ${{ANNOTATION}} -b ${{BAM_LIST}} -p ${{TMP}} " \
    #                  "-o ${{PREFIX}}\"candidate_list_from_disc.txt\" -n {1} {2}\n".format(idisc_c, ncores, s_user)
    # sbarcode_step = f"python3 ${{XTEA_PATH}}\"x_TEA_main.py\" -B -i ${{PREFIX}}\"candidate_list_from_disc.txt\" --nb 400 " \
    #                 "--ref ${{REF}} -a ${{ANNOTATION}} -b ${{BAM1}} -d ${{BARCODE_BAM}} -p ${{TMP}} " \
    #                 "-o ${{PREFIX}}\"candidate_list_barcode.txt\" -n {ncores} {s_user}\n"
    # sfilter_10x = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" -N --cr {iflt_clip} --nd {iflt_disc} -b ${{BAM_LIST}} -p ${{TMP_CNS}} " \
    #               "--fflank ${{SF_FLANK}} --flklen {iflk_len} -n {ncores} -i ${{PREFIX}}\"candidate_list_barcode.txt\" " \
    #               "-r ${{L1_CNS}} --ref ${{REF}} -a ${{ANNOTATION}} " \
    #               "-o ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" {s_user} {s_resume}\n"
    # if b_SVA==True:
    #     sfilter_10x = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" -N --sva --cr {0} --nd {1} -b ${{BAM_LIST}} -p ${{TMP_CNS}} " \
    #                   "--fflank ${{SF_FLANK}} --flklen {2} -n {3} -i ${{PREFIX}}\"candidate_list_barcode.txt\" " \
    #                   "-r ${{L1_CNS}} --ref ${{REF}} -a ${{ANNOTATION}} " \
    #                   "-o ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" {4} {5} {6} {7}\n".format(iflt_clip, iflt_disc, iflk_len,
    #                                                                                    ncores, s_purity, s_user, s_tumor, s_resume)

#     s_filter = f"python3 ${{XTEA_PATH}}\"x_TEA_main.py\" -N --cr {0} --nd {1} -b ${{BAM_LIST}} -p ${{TMP_CNS}} " \
#                "--fflank ${{SF_FLANK}} --flklen {2} -n {3} -i ${{PREFIX}}\"candidate_list_from_disc.txt\" " \
#                "-r ${{L1_CNS}} --ref ${{REF}} -a ${{ANNOTATION}} " \
#                "-o ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" {4} {5} {6} {7}\n".format(iflt_clip, iflt_disc, iflk_len, ncores, s_purity, s_user, s_tumor, s_resume)
#     if b_SVA==True:
#         s_filter = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" -N --sva --cr {0} --nd {1} -b ${{BAM_LIST}} -p ${{TMP_CNS}} " \
#                    "--fflank ${{SF_FLANK}} --flklen {2} -n {3} -i ${{PREFIX}}\"candidate_list_from_disc.txt\" " \
#                    "-r ${{L1_CNS}} --ref ${{REF}} -a ${{ANNOTATION}} " \
#                    "-o ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" {4} {5} {6} {7}\n".format(iflt_clip, iflt_disc, iflk_len,
#                                                                                     ncores, s_purity, s_user, s_tumor, s_resume)
#     if b_L1==True:
#         s_filter = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" -N --l1 --cr {0} --nd {1} -b ${{BAM_LIST}} -p ${{TMP_CNS}} " \
#                    "--fflank ${{SF_FLANK}} --flklen {2} -n {3} -i ${{PREFIX}}\"candidate_list_from_disc.txt\" " \
#                    "-r ${{L1_CNS}} --ref ${{REF}} -a ${{ANNOTATION}} " \
#                    "-o ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" {4} {5} {6} {7}\n".format(iflt_clip, iflt_disc, iflk_len,
#                                                                                     ncores, s_purity, s_user, s_tumor, s_resume)
#     sf_trsdct = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" --transduction --cr {0} --nd {1} -b ${{BAM_LIST}} -p ${{TMP_TNSD}} " \
#                 "--fflank ${{SF_FLANK}} --flklen {2} -n {3} -i ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" " \
#                 "-r ${{L1_CNS}} --ref ${{REF}} --input2 ${{PREFIX}}\"candidate_list_from_disc.txt.clip_sites_raw_disc.txt\" " \
#                 "--rtype {4} -a ${{ANNOTATION1}} {5} {6} {7} " \
#                 "-o ${{PREFIX}}\"candidate_disc_filtered_cns2.txt\"\n".format(iflt_clip, iflt_disc, iflk_len, ncores, i_rep_type, s_purity, s_tumor, s_resume)
#     sf_sibling = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" --sibling --cr {0} --nd {1} -b ${{BAM_LIST}} -p ${{TMP_TNSD}} " \
#                 "--fflank ${{SF_FLANK}} --flklen {2} -n {3} -i ${{PREFIX}}\"candidate_disc_filtered_cns2.txt\" " \
#                 "-r ${{L1_CNS}} --ref ${{REF}} --input2 ${{PREFIX}}\"candidate_list_from_disc.txt.clip_sites_raw_disc.txt\" " \
#                 "--rtype {4} -a ${{ANNOTATION1}} --blacklist ${{BLACK_LIST}} {5} {6} {7} " \
#                 "-o ${{PREFIX}}\"candidate_sibling_transduction2.txt\"\n".format(iflt_clip, iflt_disc, iflk_len, ncores,
#                                                                               i_rep_type, s_purity, s_tumor, s_resume)
#     if b_SVA == True:#here skip orphan transduction calling for SVA (as overlap with L1)
#         sf_sibling = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" --sibling --cr {0} --nd {1} -b ${{BAM_LIST}} -p ${{TMP_TNSD}} " \
#                      "--fflank \"\" --flklen {2} -n {3} -i ${{PREFIX}}\"candidate_disc_filtered_cns2.txt\" " \
#                      "-r ${{L1_CNS}} --ref ${{REF}} --input2 ${{PREFIX}}\"candidate_list_from_disc.txt.clip_sites_raw_disc.txt\" " \
#                      "--rtype {4} -a ${{ANNOTATION1}} --blacklist ${{BLACK_LIST}} {5} {6} {7} " \
#                      "-o ${{PREFIX}}\"candidate_sibling_transduction2.txt\"\n".format(iflt_clip, iflt_disc, iflk_len,
#                                                                                       ncores, i_rep_type, s_purity, s_tumor, s_resume)
#     sf_collect = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" -E --nb 500 -b ${{BAM1}} -d ${{BARCODE_BAM}} --ref ${{REF}} " \
#                  "-i ${{PREFIX}}\"candidate_disc_filtered_cns_post_filtering.txt\" -p ${{TMP}} -a ${{ANNOTATION}} -n {0} " \
#                  "--flklen {1}\n".format(ncores, iflk_len)
#     sf_asm = "python ${{XTEA_PATH}}\"x_TEA_main.py\" -A -L -p ${{TMP}} --ref ${{REF}} -n {0} " \
#              "-i ${{PREFIX}}\"candidate_disc_filtered_cns_post_filtering.txt\"\n".format(ncores)
#     sf_alg_ctg = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" -M -i ${{PREFIX}}\"candidate_disc_filtered_cns_post_filtering.txt\" " \
#                  "-r ${{REF}} -n {0} -p ${{TMP}} --ref ${{L1_CNS}} " \
#                  "-o ${{PREFIX}}\"candidate_list_asm.txt\"\n".format(ncores)
#     sf_mutation = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" -I -p ${{TMP}} -n {0} " \
#                   "-i ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" -r ${{L1_CNS}} " \
#                   "--teilen {1} -o ${{PREFIX}}\"internal_snp.vcf.gz\"\n".format(ncores, min_tei_len)
# ####
#     ####
#     sf_post_filter = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" --postF --rtype {0} -p ${{TMP_CNS}} " \
#                      "-n {1} -i ${{PREFIX}}\"candidate_disc_filtered_cns2.txt\" " \
#                      "-a ${{ANNOTATION1}} {2} " \
#                      "-o ${{PREFIX}}\"candidate_disc_filtered_cns_post_filtering.txt\"\n".format(i_rep_type, ncores, s_tumor)
# ####
#     if b_mosaic is True:
#         sf_post_filter = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" --postF --postFmosaic --rtype {0} -p ${{TMP_CNS}} " \
#                          "-n {1} -i ${{PREFIX}}\"candidate_disc_filtered_cns2.txt\" " \
#                          "-a ${{ANNOTATION1}} {2} " \
#                          "-o ${{PREFIX}}\"candidate_disc_filtered_cns_post_filtering.txt\"\n".format(i_rep_type, ncores, s_tumor)
# 
#     sf_post_filter_hc = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" --postF --rtype {0} -p ${{TMP_CNS}} -n {1} " \
#                         "-i ${{PREFIX}}\"candidate_disc_filtered_cns2.txt.high_confident\" -a ${{ANNOTATION1}} " \
#                         "--blacklist ${{BLACK_LIST}} {2} " \
#                         "-o ${{PREFIX}}\"candidate_disc_filtered_cns.txt.high_confident.post_filtering.txt\"\n" \
#         .format(i_rep_type, ncores, s_tumor)
#     if b_mosaic is True:
#         sf_post_filter_hc = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" --postF --postFmosaic --rtype {0} -p ${{TMP_CNS}} -n {1} " \
#                             "-i ${{PREFIX}}\"candidate_disc_filtered_cns2.txt.high_confident\" -a ${{ANNOTATION1}} {2} " \
#                             "-o ${{PREFIX}}\"candidate_disc_filtered_cns.txt.high_confident.post_filtering.txt\"\n" \
#             .format(i_rep_type, ncores, s_tumor)

####
    ####
    # s_case_control = "python ${{XTEA_PATH}}\"x_TEA_main.py\" --case_control -p ${{TMP_CNS}} -b ${{BAM_LIST}} " \
    #                  "-n {0} -i ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" " \
    #                  " --ref ${{REF}} --cr {1} --nd {2} " \
    #                  "-o ${{PREFIX}}\"candidate_disc_filtered_cns_somatic.txt\"\n" \
    #     .format(ncores, iflt_clip, iflt_disc)
    # s_case_control_hc = "python ${{XTEA_PATH}}\"x_TEA_main.py\" --case_control -p ${{TMP_CNS}} -b ${{BAM_LIST}} " \
    #                  "-n {0} -i ${{PREFIX}}\"candidate_disc_filtered_cns_post_filtering.txt\" " \
    #                  " --ref ${{REF}} --cr {1} --nd {2} " \
    #                  "-o ${{PREFIX}}\"candidate_disc_filtered_cns_post_filtering_somatic.txt\"\n" \
    #     .format(ncores, iflt_clip, iflt_disc)
    #
    # ####
    # s_igv="python ${{XTEA_PATH}}\"x_TEA_main.py\" --igv -p ${{PREFIX}}\"tmp/igv\" -b ${{BAM_LIST}} " \
    #                  "-i ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" " \
    #                  " --ref \"hg38\" -e 500 " \
    #                  "-o ${{PREFIX}}\"tmp/igv/igv_screenshot.txt\"\n" \
    #     .format(ncores, iflt_clip, iflt_disc)
    ######

    s_cmd = setup_step # YW 2021/05/26 changed from ""
    if iflag & 1 == 1:
        s_cmd += sclip_step
    if iflag & 2 == 2:
        s_cmd += sclip_locus_step
    if iflag & 4 == 4:
        if iflag & 1 == 1:
            s_cmd += sdisc_step
        if iflag & 2 == 2:
            s_cmd += sdisc_locus_step
    # if iflag & 4 == 4:
    #     s_cmd += sbarcode_step
    # if iflag & 8 == 8:
    #     s_cmd += sfilter_10x
    # if iflag & 16 == 16:
    #     s_cmd += s_filter
    # if iflag & 32 == 32:
    #     s_cmd += sf_collect
    # if iflag & 64 == 64:
    #     s_cmd += sf_asm
    # if iflag & 128 == 128:
    #     s_cmd += sf_alg_ctg
    # if iflag & 256 == 256:
    #     s_cmd += sf_trsdct
    # if iflag & 4096 == 4096:
    #     s_cmd+=sf_sibling
    # if iflag & 512 == 512:
    #     s_cmd += sf_post_filter
    #     s_cmd += sf_post_filter_hc
    # if (iflag & 1024 == 1024) and (iflag & 2048 != 2048):
    #     s_cmd += gnrt_calling_command_annotation_vcf(ncores, i_rep_type)
    return s_cmd

####
def gnrt_calling_command_annotation_vcf(ncores, i_rep_type):
    sf_gene = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" --gene -a ${{GENE}} -i ${{PREFIX}}\"candidate_disc_filtered_cns.txt.high_confident.post_filtering.txt\" " \
              " -n {0} -o ${{PREFIX}}\"candidate_disc_filtered_cns.txt.high_confident.post_filtering_with_gene.txt\"\n".format(ncores)

    sf_gntp_classify = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" --gntp_classify -i ${{PREFIX}}\"candidate_disc_filtered_cns.txt.high_confident.post_filtering_with_gene.txt\" " \
                       " -n {0} --model ${{XTEA_PATH}}\"genotyping/trained_model_ssc_py2_random_forest_two_category.pkl\" " \
                       " -o ${{PREFIX}}\"candidate_disc_filtered_cns.txt.high_confident.post_filtering_with_gene_gntp.txt\"\n".format(1)

    sf_cvt_gvcf = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" --gVCF -i ${{PREFIX}}\"candidate_disc_filtered_cns.txt.high_confident.post_filtering_with_gene_gntp.txt\" " \
                  " -o ${{PREFIX}} -b ${{BAM_LIST}} --ref ${{REF}} --rtype {0}\n".format(i_rep_type)
    s_cmd=""
    s_cmd += sf_gene
    s_cmd += sf_gntp_classify
    s_cmd += sf_cvt_gvcf
    return s_cmd

def gnrt_calling_command_annotation_vcf_somatic(ncores, i_rep_type):
    sf_gene = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" --gene -a ${{GENE}} -i ${{PREFIX}}\"candidate_disc_filtered_cns_high_confident_post_filtering_somatic.txt\" " \
              " -n {0} -o ${{PREFIX}}\"candidate_disc_filtered_cns.txt.high_confident.post_filtering_somatic_with_gene.txt\"\n".format(ncores)

    sf_gntp_classify = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" --gntp_classify -i ${{PREFIX}}\"candidate_disc_filtered_cns.txt.high_confident.post_filtering_somatic_with_gene.txt\" " \
                       " -n {0} --model ${{XTEA_PATH}}\"genotyping/trained_model_ssc_py2_random_forest_two_category.pkl\" " \
                       " -o ${{PREFIX}}\"candidate_disc_filtered_cns.txt.high_confident.post_filtering_somatic_with_gene_gntp.txt\"\n".format(1)

    sf_cvt_gvcf = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" --gVCF -i ${{PREFIX}}\"candidate_disc_filtered_cns.txt.high_confident.post_filtering_somatic_with_gene_gntp.txt\" " \
                  " -o ${{PREFIX}} -b ${{BAM_LIST}} --ref ${{REF}} --rtype {0}\n".format(i_rep_type)
    s_cmd=""
    s_cmd += sf_gene
    s_cmd += sf_gntp_classify
    s_cmd += sf_cvt_gvcf
    return s_cmd

####
def gnrt_calling_command_somatic_case_control(iflt_clip, iflt_disc, ncores, sf_ctr_list, iflk_len,
                                              sf_case_control_bam_list, iflag, i_rep_type):
    s_cmd = ""
    s_case_control = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" --case_control -p ${{TMP_CNS}} -b {0} " \
                     "-n {1} -r ${{L1_CNS}} -i ${{PREFIX}}\"candidate_disc_filtered_cns_post_filtering.txt\" " \
                     " --ref ${{REF}} --cr {2} --nd {3} --fflank ${{SF_FLANK}} --flklen {4} " \
                     "-o ${{PREFIX}}\"candidate_disc_filtered_cns_post_filtering_somatic.txt\"\n" \
        .format(sf_ctr_list, ncores, iflt_clip, iflt_disc, iflk_len)
    s_case_control_hc = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" --case_control --somatic_hc -p ${{TMP_CNS}} " \
                        "-i ${{PREFIX}}\"candidate_disc_filtered_cns.txt.high_confident.post_filtering.txt\" " \
                        " --input2 ${{PREFIX}}\"candidate_disc_filtered_cns_post_filtering_somatic.txt\" -n {0} " \
                        "--fflank ${{SF_FLANK}} --flklen {1} " \
                        "-o ${{PREFIX}}\"candidate_disc_filtered_cns_high_confident_post_filtering_somatic.txt\"\n".format(ncores, iflk_len)
    ####
    s_igv="python3 ${{XTEA_PATH}}\"x_TEA_main.py\" --igv --single_sample -p ${{PREFIX}}\"tmp/igv\" -b {0} " \
                     "-i ${{PREFIX}}\"candidate_disc_filtered_cns2.txt\" " \
                     " --ref \"hg38\" -e 500 " \
                     "-o ${{PREFIX}}\"tmp/igv/igv_screenshot.txt\"\n".format(sf_case_control_bam_list)
    s_igv_hc = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" --igv --single_sample -p ${{PREFIX}}\"tmp/igv\" -b {0} " \
            "-i ${{PREFIX}}\"candidate_disc_filtered_cns_high_confident_post_filtering_somatic.txt\" " \
            " --ref \"hg38\" -e 500 " \
            "-o ${{PREFIX}}\"tmp/igv/igv_screenshot_somatic.txt\"\n".format(sf_case_control_bam_list)

    ####BamSnap
    s_bamsnap = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" --igv --bamsnap --single_sample -p ${{PREFIX}}\"tmp/igv\" -b {0} " \
            "-i ${{PREFIX}}\"candidate_disc_filtered_cns2.txt\" " \
            " --ref ${{REF}} -e 1000 " \
            "-o ${{PREFIX}}\"tmp/igv/bamsnap_screenshot.txt\"\n".format(sf_case_control_bam_list)
    s_bamsnap_hc = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" --igv --bamsnap --single_sample -p ${{PREFIX}}\"tmp/igv\" -b {0} " \
               "-i ${{PREFIX}}\"candidate_disc_filtered_cns_high_confident_post_filtering_somatic.txt\" " \
               " --ref ${{REF}} -e 1000 " \
               "-o ${{PREFIX}}\"tmp/igv/bamsnap_screenshot_somatic.txt\"\n".format(sf_case_control_bam_list)
    ####
    s_cmd += s_case_control
    s_cmd += s_case_control_hc
    s_cmd += s_igv
    s_cmd += s_igv_hc
    s_cmd += s_bamsnap
    s_cmd += s_bamsnap_hc

    if iflag & 1024 == 1024:
        s_cmd += gnrt_calling_command_annotation_vcf_somatic(ncores, i_rep_type)
    return s_cmd
####
####
def gnrt_calling_command_non_RNA(iclip_c, iclip_rp, idisc_c, iflt_clip, iflt_disc, ncores, iflk_len, min_tei_len, iflag,
                                 b_mosaic, b_user_par, s_cfolder):
    s_user = ""
    if b_user_par == True:
        s_user = "--user"
    sclip_step = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" -C --dna -i ${{BAM_LIST}} --lc {0} --rc {1} --cr {2}  " \
                 "-r ${{L1_COPY_WITH_FLANK}}  -a ${{ANNOTATION}} --cns ${{L1_CNS}} --ref ${{REF}} -p ${{TMP}} " \
                 "-o ${{PREFIX}}\"candidate_list_from_clip.txt\"  -n {3} --cp {4} {5}\n".format(iclip_c, iclip_c, iclip_rp,
                                                                                            ncores, s_cfolder, s_user)
    sdisc_step = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\"  -D --dna -i ${{PREFIX}}\"candidate_list_from_clip.txt\" --nd {0} " \
                 "--ref ${{REF}} -a ${{ANNOTATION}} -b ${{BAM_LIST}} -p ${{TMP}} " \
                 "-o ${{PREFIX}}\"candidate_list_from_disc.txt\" -n {1} {2}\n".format(idisc_c, ncores, s_user)
    if b_mosaic == True:
        sdisc_step = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\"  -D --somatic --dna -i ${{PREFIX}}\"candidate_list_from_clip.txt\" --nd {0} " \
                     "--ref ${{REF}} -a ${{ANNOTATION}} -b ${{BAM_LIST}} -p ${{TMP}} " \
                     "-o ${{PREFIX}}\"candidate_list_from_disc.txt\" -n {1} {2}\n".format(idisc_c, ncores, s_user)
    sbarcode_step = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" -B --dna -i ${{PREFIX}}\"candidate_list_from_disc.txt\" --nb 400 " \
                    "--ref ${{REF}} -a ${{ANNOTATION}} -b ${{BAM1}} -d ${{BARCODE_BAM}} -p ${{TMP}} " \
                    "-o ${{PREFIX}}\"candidate_list_barcode.txt\" -n {0} {1}\n".format(ncores, s_user)
    sfilter_10x = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" -N --dna --cr {0} --nd {1} -b ${{BAM_LIST}} -p ${{TMP_CNS}} " \
                  "--fflank ${{SF_FLANK}} --flklen {2} -n {3} -i ${{PREFIX}}\"candidate_list_barcode.txt\" " \
                  "-r ${{L1_CNS}} --ref ${{REF}} -a ${{ANNOTATION}} " \
                  "-o ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" {4}\n".format(iflt_clip, iflt_disc, iflk_len, ncores, s_user)
    s_filter = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" -N --dna --cr {0} --nd {1} -b ${{BAM_LIST}} -p ${{TMP_CNS}} " \
               "--fflank ${{SF_FLANK}} --flklen {2} -n {3} -i ${{PREFIX}}\"candidate_list_from_disc.txt\" " \
               "-r ${{L1_CNS}} --ref ${{REF}} -a ${{ANNOTATION}} " \
               "-o ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" {4}\n".format(iflt_clip, iflt_disc, iflk_len, ncores, s_user)
    sf_collect = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" -E --dna --nb 500 -b ${{BAM1}} -d ${{BARCODE_BAM}} --ref ${{REF}} " \
                 "-i ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" -p ${{TMP}} -a ${{ANNOTATION}} -n {0} " \
                 "--flklen {1}\n".format(ncores, iflk_len)
    sf_asm = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" -A -L --dna -p ${{TMP}} --ref ${{REF}} -n {0} " \
             "-i ${{PREFIX}}\"candidate_disc_filtered_cns.txt\"\n".format(ncores)
    sf_alg_ctg = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" -M --dna -i ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" " \
                 "-r ${{REF}} -n {0} -p ${{TMP}} --ref ${{L1_CNS}} " \
                 "-o ${{PREFIX}}\"candidate_list_asm.txt\"\n".format(ncores)
    sf_mutation = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" -I --dna -p ${{TMP}} -n {0} " \
                  "-i ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" -r ${{L1_CNS}} " \
                  "--teilen {1} -o ${{PREFIX}}\"internal_snp.vcf.gz\"\n".format(ncores, min_tei_len)
    sf_gene="python3 ${{XTEA_PATH}}\"x_TEA_main.py\" --gene -a ${{GENE}} -i ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" " \
            " -n {0} -o ${{PREFIX}}\"candidate_disc_filtered_cns_with_gene.txt\"\n".format(ncores)
    #sf_clean_tmp = "find ${TMP} -type f -name \'*tmp*\' -delete\n"
    # sf_clean_sam="find ${TMP} -type f -name \'*.sam\' -delete\n"
    # sf_clean_fa="find ${TMP} -type f -name \'*.fa\' -delete\n"
    # sf_clean_fq = "find ${TMP} -type f -name \'*.fq\' -delete\n"

    ####
    ####
    s_cmd = ""
    if iflag & 1 == 1:
        s_cmd += sclip_step
    if iflag & 2 == 2:
        s_cmd += sdisc_step
    if iflag & 4 == 4:
        s_cmd += sbarcode_step
    if iflag & 8 == 8:
        s_cmd += sfilter_10x
    if iflag & 16 == 16:
        s_cmd += s_filter
    if iflag & 32 == 32:
        s_cmd += sf_collect
    if iflag & 64 == 64:
        s_cmd += sf_asm
    if iflag & 128 == 128:
        s_cmd += sf_alg_ctg
    if iflag & 256 == 256:
        s_cmd += sf_mutation
    if iflag & 512 == 512:
        s_cmd += sf_gene

    # s_cmd+=sf_clean_tmp
    # s_cmd +=sf_clean_sam
    # s_cmd +=sf_clean_fa
    # s_cmd +=sf_clean_fq
    return s_cmd

# grnt calling steps
def gnrt_calling_command_MT(iclip_c, iclip_rp, idisc_c, iflt_clip, iflt_disc, ncores, iflk_len, min_tei_len, iflag,
                            b_mosaic, b_user_par, s_cfolder):
    s_user = ""
    if b_user_par == True:
        s_user = "--user"
    sclip_step = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" -C --dna --mit -i ${{BAM_LIST}} --lc {0} --rc {1} --cr {2}  " \
                 "-r ${{L1_COPY_WITH_FLANK}}  -a ${{ANNOTATION}} --cns ${{L1_CNS}} --ref ${{REF}} -p ${{TMP}} " \
                 "-o ${{PREFIX}}\"candidate_list_from_clip.txt\"  -n {3} --cp {4} {5}\n".format(iclip_c, iclip_c, iclip_rp,
                                                                                            ncores, s_cfolder, s_user)
    sdisc_step = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\"  -D --dna --mit -i ${{PREFIX}}\"candidate_list_from_clip.txt\" --nd {0} " \
                 "--ref ${{REF}} -a ${{ANNOTATION}} -b ${{BAM_LIST}} -p ${{TMP}} " \
                 "-o ${{PREFIX}}\"candidate_list_from_disc.txt\" -n {1} {2}\n".format(idisc_c, ncores, s_user)
    if b_mosaic == True:
        sdisc_step = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\"  -D --somatic --dna --mit -i ${{PREFIX}}\"candidate_list_from_clip.txt\" --nd {0} " \
                     "--ref ${{REF}} -a ${{ANNOTATION}} -b ${{BAM_LIST}} -p ${{TMP}} " \
                     "-o ${{PREFIX}}\"candidate_list_from_disc.txt\" -n {1} {2}\n".format(idisc_c, ncores, s_user)
    sbarcode_step = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" -B --dna --mit -i ${{PREFIX}}\"candidate_list_from_disc.txt\" --nb 400 " \
                    "--ref ${{REF}} -a ${{ANNOTATION}} -b ${{BAM1}} -d ${{BARCODE_BAM}} -p ${{TMP}} " \
                    "-o ${{PREFIX}}\"candidate_list_barcode.txt\" -n {0} {1}\n".format(ncores, s_user)
    sfilter_10x = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" -N --dna --mit --cr {0} --nd {1} -b ${{BAM_LIST}} -p ${{TMP_CNS}} " \
                  "--fflank ${{SF_FLANK}} --flklen {2} -n {3} -i ${{PREFIX}}\"candidate_list_barcode.txt\" " \
                  "-r ${{L1_CNS}} --ref ${{REF}} -a ${{ANNOTATION}} " \
                  "-o ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" {4}\n".format(iflt_clip, iflt_disc, iflk_len, ncores, s_user)
    s_filter = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" -N --dna --mit --cr {0} --nd {1} -b ${{BAM_LIST}} -p ${{TMP_CNS}} " \
               "--fflank ${{SF_FLANK}} --flklen {2} -n {3} -i ${{PREFIX}}\"candidate_list_from_disc.txt\" " \
               "-r ${{L1_CNS}} --ref ${{REF}} -a ${{ANNOTATION}} " \
               "-o ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" {4}\n".format(iflt_clip, iflt_disc, iflk_len, ncores, s_user)
    sf_collect = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" -E --dna --mit --nb 500 -b ${{BAM1}} -d ${{BARCODE_BAM}} --ref ${{REF}} " \
                 "-i ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" -p ${{TMP}} -a ${{ANNOTATION}} -n {0} " \
                 "--flklen {1}\n".format(ncores, iflk_len)
    sf_asm = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" -A -L --dna --mit -p ${{TMP}} --ref ${{REF}} -n {0} " \
             "-i ${{PREFIX}}\"candidate_disc_filtered_cns.txt\"\n".format(ncores)
    sf_alg_ctg = "python3 ${{XTEA_PATH}}\"x_TEA_main.py\" -M --dna --mit -i ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" " \
                 "-r ${{REF}} -n {0} -p ${{TMP}} --ref ${{L1_CNS}} " \
                 "-o ${{PREFIX}}\"candidate_list_asm.txt\"\n".format(ncores)
    sf_mutation = "python ${{XTEA_PATH}}\"x_TEA_main.py\" -I --dna --mit -p ${{TMP}} -n {0} " \
                  "-i ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" -r ${{L1_CNS}} " \
                  "--teilen {1} -o ${{PREFIX}}\"internal_snp.vcf.gz\"\n".format(ncores, min_tei_len)
    sf_gene="python ${{XTEA_PATH}}\"x_TEA_main.py\" --gene --dna --mit -a ${{GENE}} -i ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" " \
            " -n {0} -o ${{PREFIX}}\"candidate_disc_filtered_cns_with_gene.txt\"\n".format(ncores)
    #sf_clean_tmp = "find ${TMP} -type f -name \'*tmp*\' -delete\n"
    # sf_clean_sam="find ${TMP} -type f -name \'*.sam\' -delete\n"
    # sf_clean_fa="find ${TMP} -type f -name \'*.fa\' -delete\n"
    # sf_clean_fq = "find ${TMP} -type f -name \'*.fq\' -delete\n"

    ####
    ####
    s_cmd = ""
    if iflag & 1 == 1:
        s_cmd += sclip_step
    if iflag & 2 == 2:
        s_cmd += sdisc_step
    if iflag & 4 == 4:
        s_cmd += sbarcode_step
    if iflag & 8 == 8:
        s_cmd += sfilter_10x
    if iflag & 16 == 16:
        s_cmd += s_filter
    if iflag & 32 == 32:
        s_cmd += sf_collect
    if iflag & 64 == 64:
        s_cmd += sf_asm
    if iflag & 128 == 128:
        s_cmd += sf_alg_ctg
    if iflag & 256 == 256:
        s_cmd += sf_mutation
    if iflag & 512 == 512:
        s_cmd += sf_gene

    # s_cmd+=sf_clean_tmp
    # s_cmd +=sf_clean_sam
    # s_cmd +=sf_clean_fa
    # s_cmd +=sf_clean_fq
    return s_cmd


####gnrt the whole pipeline
# YW 2021/05/23 removed rep_type, use the sf_root_folder as it is, not sf_working_folder, commented out code to generate cns, tsdct, igv folders
def gnrt_pipelines(s_head, s_libs, s_calling_cmd, sf_id, sf_bams, sf_bams_10X, sf_root_folder):
    l_sbatch_files=[]
    # sf_working_folder=sf_root_folder+rep_type
    if sf_root_folder[-1] != "/":
        sf_root_folder += "/"

    m_id = {}
    with open(sf_id) as fin_id:
        for line in fin_id:
            sid = line.rstrip()
            m_id[sid] = 1

            sf_folder = sf_root_folder
            if os.path.exists(sf_folder)==False:
                cmd = "mkdir {0}".format(sf_folder)
                run_cmd(cmd)
            # create the temporary folders
            sf_tmp=sf_folder + "/tmp"
            if os.path.exists(sf_tmp)==False:
                cmd = "mkdir {0}".format(sf_tmp)
                run_cmd(cmd)
            sf_tmp_clip=sf_folder + "/tmp/clip"
            if os.path.exists(sf_tmp_clip)==False:
                cmd = "mkdir {0}".format(sf_tmp_clip)
                run_cmd(cmd)
            # sf_tmp_cns=sf_folder + "/tmp/cns"
            # if os.path.exists(sf_tmp_cns)==False:
            #     cmd = "mkdir {0}".format(sf_tmp_cns)
            #     run_cmd(cmd)
            # sf_tmp_tsdct = sf_folder + "/tmp/transduction"
            # if os.path.exists(sf_tmp_tsdct) == False:
            #     cmd = "mkdir {0}".format(sf_tmp_tsdct)
            #     run_cmd(cmd)
            # sf_tmp_igv = sf_folder + "/tmp/igv"
            # if os.path.exists(sf_tmp_igv) == False:
            #     cmd = "mkdir {0}".format(sf_tmp_igv)
            #     run_cmd(cmd)
    m_bams = {}
    if sf_bams != "null":
        with open(sf_bams) as fin_bams:
            for line in fin_bams:
                fields = line.split()
                sid = fields[0]
                s_bam = fields[1]
                if sid not in m_bams:
                    m_bams[sid]=[]
                m_bams[sid].append(s_bam)

    m_bams_10X = {}
    if sf_bams_10X != "null":
        with open(sf_bams_10X) as fin_bams_10X:
            for line in fin_bams_10X:
                fields = line.split()
                sid = fields[0]
                if sid not in m_id:
                    continue
                s_bam = fields[1]
                s_barcode_bam = fields[2]
                m_bams_10X[sid] = s_bam

                # soft-link the bams
                sf_10X_bam = sf_working_folder + "10X_phased_possorted_bam.bam"
                if os.path.isfile(sf_10X_bam) == False:
                    cmd = f"ln -s {s_bam} {sf_10X_bam}"
                    run_cmd(cmd)

                sf_10X_barcode_bam = sf_working_folder + "10X_barcode_indexed.sorted.bam"
                if os.path.isfile(sf_10X_barcode_bam) == False:
                    cmd = f"ln -s {s_barcode_bam} {sf_10X_barcode_bam}"
                    run_cmd(cmd)
                # soft-link the bai
                sf_10X_bai = sf_working_folder + "10X_phased_possorted_bam.bam.bai"
                if os.path.isfile(sf_10X_bai) == False:
                    cmd = f"ln -s {s_bam + '.bai'} {sf_10X_bai}"
                    run_cmd(cmd)

                sf_10X_barcode_bai = sf_working_folder + "10X_barcode_indexed.sorted.bam.bai"
                if os.path.isfile(sf_10X_barcode_bai) == False:
                    cmd = f"ln -s {s_barcode_bam + '.bai'} {sf_10X_barcode_bai}"
                    run_cmd(cmd)
    
####
    for sid in m_id:
        sf_folder = sf_root_folder
        if os.path.exists(sf_folder) == False:
            continue

        ####gnrt the bam list file
        sf_bam_list = sf_folder + "bam_list.txt"
        with open(sf_bam_list, "w") as fout_bam_list:
            if sid in m_bams:
                for sf_tmp_bam in m_bams[sid]:
                    fout_bam_list.write(sf_tmp_bam + "\t"+ ILLUMINA + "\n")
            if sid in m_bams_10X:
                fout_bam_list.write(m_bams_10X[sid] + "\t"+ X10 + "\n")
        ####gnrt the pipeline file
        sf_out_sh = sf_folder + "run_xTEA_pipeline.sh"
        with open(sf_out_sh, "w") as fout_sh:  ###gnrt the pipeline file
            fout_sh.write(s_head)
            s_prefix = "PREFIX={0}\n".format(sf_folder)
            fout_sh.write(s_prefix)
            fout_sh.write("############\n")
            fout_sh.write("############\n")
            fout_sh.write(s_libs)
            fout_sh.write("############\n")
            fout_sh.write("############\n")
            fout_sh.write(s_calling_cmd)
        l_sbatch_files.append(sf_out_sh)

    return l_sbatch_files

# YW 2021/05/20 took out s_tmp_cns, sf_anno1, sf_flank; substitute sf_anno, sf_copy_with_flank, sf_cns with l_sf_anno, l_sf_copy_with_flank, l_sf_cns respectively
# YW 2021/08/04 added sf_ctrl_bed
def write_to_config(l_sf_anno, sf_ref, sf_gene, sf_black_list, l_sf_copy_with_flank, sf_flank, l_sf_cns,
                    sf_xtea, s_bl, s_bam1, s_bc_bam, s_tmp, s_tmp_clip, sf_ctrl_bed, sf_config):
    with open(sf_config, "w") as fout:
        for sf_anno in l_sf_anno:
            fout.write(sf_anno)
        fout.write(sf_ref)
        fout.write(sf_gene)
        fout.write(sf_black_list)
        for sf_copy_with_flank in l_sf_copy_with_flank:
            fout.write(sf_copy_with_flank)
        fout.write(sf_flank)
        for sf_cns in l_sf_cns:
            fout.write(sf_cns)
        fout.write(sf_xtea)
        fout.write(s_bl)
        fout.write(s_bam1)
        fout.write(s_bc_bam)
        fout.write(s_tmp)
        fout.write(s_tmp_clip)
        fout.write(sf_ctrl_bed)
        # fout.write(s_tmp_cns)
#
#generate library configuration files
# YW 2021/05/19 enabled adding multiple library configuration files to the same command (separated by repeat type/MT)
# YW 2021/08/04 added sf_ctrl_bed to enable control bam file coordinate lifting
def gnrt_lib_config(l_rep_type, sf_folder_rep, sf_ref, sf_gene, sf_black_list, sf_folder_xtea, sf_config_prefix, sf_ctrl_bed):
    if sf_folder_rep[-1] != "/":
        sf_folder_rep += "/"
    if sf_folder_xtea[-1] != "/":
        sf_folder_xtea += "/"
    if sf_config_prefix[-1] != "/":
        sf_config_prefix += "/"

    s_bl = "BAM_LIST ${PREFIX}\"bam_list.txt\"\n"
    s_bam1 = "BAM1 ${PREFIX}\"10X_phased_possorted_bam.bam\"\n"
    s_bc_bam = "BARCODE_BAM ${PREFIX}\"10X_barcode_indexed.sorted.bam\"\n"
    s_tmp = "TMP ${PREFIX}\"tmp/\"\n"
    s_tmp_clip = "TMP_CLIP ${PREFIX}\"tmp/clip/\"\n"
    # s_tmp_cns = "TMP_CNS ${PREFIX}\"tmp/cns/\"\n"
    # s_tmp_tsdct = "TMP_TNSD ${PREFIX}\"tmp/transduction/\"\n"
    # s_tmp_cns += s_tmp_tsdct
    sf_ref = "REF " + sf_ref + "\n"
    sf_xtea = "XTEA_PATH " + sf_folder_xtea + "\n"
    sf_gene_anno="GENE " + sf_gene + "\n"
    sf_black_list = "BLACK_LIST " + sf_black_list + "\n"
    sf_flank = "SF_FLANK null\n"
    sf_ctrl_bed = "CTRL_BED " + sf_ctrl_bed + "\n"
    l_sf_anno = []
    # l_sf_anno1 = []
    l_sf_copy_with_flank = []
    l_sf_cns = []
    for s_rep_type in l_rep_type:
        if s_rep_type=="Alu":# for Alu
            # sf_config_Alu = sf_config_prefix + "Alu.config"
            l_sf_anno.append("ALU_ANNOTATION " + sf_folder_rep + "Alu/hg38/hg38_Alu.out\n")
            # l_sf_anno1.append("ALU_ANNOTATION " + sf_folder_rep + "Alu/hg38/hg38_Alu.out\n")
            l_sf_copy_with_flank.append("ALU_COPY_WITH_FLANK " + sf_folder_rep + "Alu/hg38/hg38_AluJabc_copies_with_flank.fa\n")
            # l_sf_flank.append("SF_FLANK_Alu null\n")
            l_sf_cns.append("ALU_CNS " + sf_folder_rep + "consensus/ALU.fa\n")
            # write_to_config(sf_anno, sf_anno1, sf_ref, sf_gene_anno, sf_black_list, sf_copy_with_flank, sf_flank,
                            # sf_cns, sf_xtea, s_bl, s_bam1, s_bc_bam, s_tmp, s_tmp_clip, s_tmp_cns, sf_config_Alu)
        elif s_rep_type=="L1":# for L1
            # sf_config_L1 = sf_config_prefix + "L1.config"
            l_sf_anno.append("L1_ANNOTATION " + sf_folder_rep + "L1/hg38/hg38_L1_larger_500_with_all_L1HS.out\n")
            # l_sf_anno1.append("L1_ANNOTATION1 " + sf_folder_rep + "L1/hg38/hg38_L1.fa.out\n")
            # YW: 2020/04/22 will include all the recently active L1 families, so change the file name
            #sf_copy_with_flank = "L1_COPY_WITH_FLANK " + sf_folder_rep + "LINE/hg38/hg38_L1HS_copies_larger_5K_with_flank.fa\n"
            l_sf_copy_with_flank.append("L1_COPY_WITH_FLANK " + sf_folder_rep + "L1/hg38/hg38_L1_copies_larger_5K_with_flank.fa\n")
            # changing this on 2020/04/22, for ancient genomes, we don't care about transduction
            #sf_flank = "SF_FLANK " + sf_folder_rep + "LINE/hg38/hg38_FL_L1_flanks_3k.fa\n"
            # l_sf_flank.append("SF_FLANK_L1 null\n")
            l_sf_cns.append("L1_CNS " + sf_folder_rep + "consensus/LINE1.fa\n")
            # write_to_config(sf_anno, sf_anno1, sf_ref, sf_gene_anno, sf_black_list, sf_copy_with_flank, sf_flank,
                            # sf_cns, sf_xtea, s_bl, s_bam1, s_bc_bam, s_tmp, s_tmp_clip, sf_config_L1)
        elif s_rep_type=="SVA":# for SVA
            # sf_config_SVA = sf_config_prefix + "SVA.config"
            l_sf_anno.append("SVA_ANNOTATION " + sf_folder_rep + "SVA/hg38/hg38_SVA.out\n")
            # l_sf_anno1.append("SVA_ANNOTATION " + sf_folder_rep + "SVA/hg38/hg38_SVA.out\n")
            l_sf_copy_with_flank.append("SVA_COPY_WITH_FLANK " + sf_folder_rep + "SVA/hg38/hg38_SVA_copies_with_flank.fa\n")
            # YW: changing this on 2020/04/22, for ancient genomes, we don't care about transduction
            #sf_flank = "SF_FLANK " + sf_folder_rep + "SVA/hg38/hg38_FL_SVA_flanks_3k.fa\n"
            # l_sf_flank.append("SF_FLANK_SVA null\n")
            l_sf_cns.append("SVA_CNS " + sf_folder_rep + "consensus/SVA.fa\n")
            # write_to_config(sf_anno, sf_anno1, sf_ref, sf_gene_anno, sf_black_list, sf_copy_with_flank, sf_flank,
                            # sf_cns, sf_xtea, s_bl, s_bam1, s_bc_bam, s_tmp, s_tmp_clip, s_tmp_cns, sf_config_SVA)
        elif s_rep_type=="HERV":####HERV
            # sf_config_HERV = sf_config_prefix + "HERV.config"
            l_sf_anno.append("HERV_ANNOTATION " + sf_folder_rep + "HERV/hg38/hg38_HERV.out\n")
            # l_sf_anno1.append("ANNOTATION1_HERV " + sf_folder_rep + "HERV/hg38/hg38_HERV.out\n")
            l_sf_copy_with_flank.append("HERV_COPY_WITH_FLANK " + sf_folder_rep + "HERV/hg38/hg38_HERV_copies_with_flank.fa\n")
            # l_sf_flank.append("SF_FLANK_HERV null\n")
            l_sf_cns.append("HERV_CNS " + sf_folder_rep + "consensus/HERV.fa\n")
            # write_to_config(sf_anno, sf_ref, sf_anno1, sf_gene_anno, sf_black_list, sf_copy_with_flank, sf_flank,
                            # sf_cns, sf_xtea, s_bl, s_bam1, s_bc_bam, s_tmp, s_tmp_clip, s_tmp_cns, sf_config_HERV)
        elif s_rep_type=="MSTA":####MSTA
            # sf_config_MSTA = sf_config_prefix + "MSTA.config"
            l_sf_anno.append("MSTA_ANNOTATION " + sf_folder_rep + "MSTA/hg38/hg38_MSTA.out\n")
            # l_sf_anno1.append("ANNOTATION1_MSTA " + sf_folder_rep + "MSTA/hg38/hg38_MSTA.out\n")
            l_sf_copy_with_flank.append("MSTA_COPY_WITH_FLANK " + sf_folder_rep + "MSTA/hg38/hg38_MSTA_copies_with_flank.fa\n")
            # l_sf_flank.append("SF_FLANK_MSTA null\n")
            l_sf_cns.append("MSTA_CNS " + sf_folder_rep + "consensus/MSTA.fa\n")
            # write_to_config(sf_anno, sf_anno1, sf_ref, sf_gene_anno, sf_black_list, sf_copy_with_flank, sf_flank,
                            # sf_cns, sf_xtea, s_bl, s_bam1, s_bc_bam, s_tmp, s_tmp_clip, s_tmp_cns, sf_config_MSTA)

        elif s_rep_type=="Mitochondrion":####for Mitochondrion
            # sf_config_Mitochondrion = sf_config_prefix + "Mitochondrion.config"
            l_sf_anno.append("MT_ANNOTATION " + sf_folder_rep + "Mitochondrion/hg38/hg38_numtS.out\n")
            # l_sf_anno1.append("ANNOTATION1_MT " + sf_folder_rep + "Mitochondrion/hg38/hg38_numtS.out\n")
            l_sf_copy_with_flank.append("MT_COPY_WITH_FLANK " + sf_folder_rep +\
                                 "Mitochondrion/hg38/hg38_mitochondrion_copies_with_flank.fa\n")
            # l_sf_flank.append("SF_FLANK_MT null\n")
            l_sf_cns.append("MT_CNS " + sf_folder_rep + "consensus/mitochondrion.fa\n")
            # write_to_config(sf_anno, sf_anno1, sf_ref, sf_gene_anno, sf_black_list, sf_copy_with_flank, sf_flank,
                            # sf_cns, sf_xtea, s_bl, s_bam1, s_bc_bam, s_tmp, s_tmp_clip, s_tmp_cns, sf_config_Mitochondrion)
    # YW 2021/05/20
    sf_config = sf_config_prefix + "config"
    write_to_config(l_sf_anno, sf_ref, sf_gene_anno, sf_black_list, l_sf_copy_with_flank, sf_flank,
                    l_sf_cns, sf_xtea, s_bl, s_bam1, s_bc_bam, s_tmp, s_tmp_clip, sf_ctrl_bed, sf_config) # YW 2021/08/04 added sf_ctrl_bed

def cp_file(sf_from, sf_to):
    cmd = "cp {0} {1}".format(sf_from, sf_to)
    if os.path.isfile(sf_from)==False:
        return
    run_cmd(cmd)

def run_cmd(cmd):
    print(cmd)
    Popen(cmd, shell=True, stdout=PIPE).communicate()

def get_sample_id(sf_bam):
    fname = ntpath.basename(sf_bam)
    fname_fields = fname.split(".")
    if fname_fields[-1] != "bam" and fname_fields[-1] != "cram":
        print("Alignment is not end with .bam")
        return None
    sample_id = ".".join(fname_fields[:-1])
    return sample_id


####gnrt the running shell
# YW 2021/05/20 modified to generate calling command using one config file, not by different rep_type, removed b_mosaic, b_tumor, f_purity
# YW 2021/08/04 added b_ctrl=False, sf_ctrl_bed="null" to enable control bam file coordinate lifting
def gnrt_running_shell(sf_ids, sf_bams, sf_10X_bams, l_rep_type, b_user_par, b_force, s_wfolder,
                       sf_folder_rep, sf_ref, sf_gene, sf_black_list, sf_folder_xtea, spartition, stime, smemory,
                       ncores, email_user, sf_submit_sh, sf_case_control_bam_list="null", b_lsf=False, b_slurm=True, b_resume=False,
                       b_ctrl=False, sf_ctrl_bed="null", sf_locus_file="null"):
    if s_wfolder[-1] != "/":
        s_wfolder += "/"
    scmd = "mkdir -p {0}".format(s_wfolder)
    Popen(scmd, shell=True, stdout=PIPE).communicate()

    m_id = {}
    with open(sf_ids) as fin_id:
        for line in fin_id:
            sid = line.rstrip()
            m_id[sid] = 1
            sf_folder = s_wfolder + sid  # first create folder
            if os.path.exists(sf_folder) == True:
                continue
            cmd = "mkdir -p {0}".format(sf_folder)
            Popen(cmd, shell=True, stdout=PIPE).communicate()
            sf_pub_clip = sf_folder+ "/" + PUB_CLIP + "/"
            cmd = "mkdir -p {0}".format(sf_pub_clip)
            Popen(cmd, shell=True, stdout=PIPE).communicate()
            # YW 2021/10/20 added the following for clip_locus_specific
            sf_pub_clip_locus = sf_folder + "/" + PUB_CLIP_LOCUS + "/"
            cmd = "mkdir -p {0}".format(sf_pub_clip_locus)
            Popen(cmd, shell=True, stdout=PIPE).communicate()

    ####gnrt the sample, bam, x10 bam files
    # YW 2021/05/19 took out the l_rep_type argument
    split_bam_list(m_id, sf_bams, sf_10X_bams, s_wfolder, sf_case_control_bam_list)
    
    # YW 2021/08/04 added to enable control bam file coordinate lifting, generate the ref_bed file
    m_ctrls = {}
    split_ctrl_bed_list(m_id, sf_ctrl_bed, s_wfolder, m_ctrls)
    
    l_sh=[]
    for sid_tmp in m_id:
        sf_sample_folder=s_wfolder + sid_tmp + "/"
        sf_pub_clip = sf_sample_folder + PUB_CLIP + "/"
        sf_pub_clip_locus = sf_sample_folder + PUB_CLIP_LOCUS + "/"
        # YW 2021/08/04 added m_ctrls to enable control bam file coordinate lifting
        if sid_tmp in m_ctrls:
            gnrt_lib_config(l_rep_type, sf_folder_rep, sf_ref, sf_gene, sf_black_list, sf_folder_xtea, sf_sample_folder, m_ctrls[sid_tmp])
        else:
            gnrt_lib_config(l_rep_type, sf_folder_rep, sf_ref, sf_gene, sf_black_list, sf_folder_xtea, sf_sample_folder, "null")
        # for rep_type in l_rep_type:
        #     i_rep_type=get_flag_by_rep_type(rep_type)
        #     sf_config = sf_sample_folder+rep_type+".config"
        # 
        #     #s_wfolder_rep=s_wfolder+rep_type
        #     s_wfolder_rep = sf_sample_folder + rep_type
        #     if os.path.exists(s_wfolder_rep)==False:
        #         cmd="mkdir {0}".format(s_wfolder_rep)
        #         run_cmd(cmd)
        # 
        #     sf_rep_sample_id = s_wfolder_rep + "/sample_id.txt"
        #     sf_rep_bam = s_wfolder_rep + "/bam_list1.txt"
        #     if os.path.isfile(sf_rep_bam)==False:
        #         sf_rep_bam="null"
        #     sf_rep_x10_bam = s_wfolder_rep + "/x10_bam_list1.txt"
        #     if os.path.exists(sf_rep_x10_bam)==False:
        #         sf_rep_x10_bam="null"
        # 
        #     s_head = ""
        #     if b_lsf==True:
        #         s_head = gnrt_script_head_lsf(spartition, ncores, stime, smemory, sid_tmp, rep_type)
        #     elif b_slurm==True:
        #         s_head = gnrt_script_head(spartition, ncores, stime, smemory, sid_tmp, rep_type)
        # 
        #     l_libs = load_par_config(sf_config)
        #     s_libs = gnrt_parameters(l_libs)
        #     ##
        #     iclip_c = args.nclip
        #     iclip_rp = args.cliprep
        #     idisc_c = args.ndisc
        #     iflt_clip = args.nfilterclip
        #     iflt_disc = args.nfilterdisc
        #     iflk_len = args.flklen
        #     itei_len = args.teilen
        #     iflag = args.flag
        # 
        #     s_calling_cmd = gnrt_calling_command(iclip_c, iclip_rp, idisc_c, iflt_clip, iflt_disc, ncores, iflk_len,
        #                                          itei_len, iflag, b_mosaic, b_user_par, b_force, b_tumor, b_resume, f_purity, i_rep_type, sf_pub_clip)
        #     if rep_type is REP_TYPE_SVA:
        #         s_calling_cmd = gnrt_calling_command(iclip_c, iclip_rp, idisc_c, iflt_clip, iflt_disc, ncores, iflk_len,
        #                                              itei_len, iflag, b_mosaic, b_user_par, b_force, b_tumor, b_resume, f_purity, i_rep_type,sf_pub_clip, True, False)
        #     # YW 2020/08/16 added the REP_TYPE_L1 statement
        #     if rep_type is REP_TYPE_L1:
        #         s_calling_cmd = gnrt_calling_command(iclip_c, iclip_rp, idisc_c, iflt_clip, iflt_disc, ncores, iflk_len,
        #                                              itei_len, iflag, b_mosaic, b_user_par, b_force, b_tumor, b_resume, f_purity, i_rep_type,sf_pub_clip, False, True)
        #     if rep_type is REP_TYPE_HERV or rep_type is REP_TYPE_MSTA:
        #         s_calling_cmd = gnrt_calling_command_non_RNA(iclip_c, iclip_rp, idisc_c, iflt_clip, iflt_disc, ncores, iflk_len,
        #                                              itei_len, iflag, b_mosaic, b_user_par, sf_pub_clip)
        #     if rep_type is REP_TYPE_MIT:
        #         s_calling_cmd = gnrt_calling_command_MT(iclip_c, iclip_rp, idisc_c, iflt_clip, iflt_disc, ncores, iflk_len,
        #                                              itei_len, iflag, b_mosaic, b_user_par, sf_pub_clip)
        #     
        #     if sf_case_control_bam_list != "null":
        #         sf_ctr_list=s_wfolder_rep + "/bam_list2.txt"
        #         sf_case_ctrl_list=s_wfolder_rep + "/bam_list3.txt"
        #         sf_tmp_cmd=gnrt_calling_command_somatic_case_control(iflt_clip, iflt_disc, ncores, sf_ctr_list,
        #                                                              iflk_len, sf_case_ctrl_list, iflag, i_rep_type)
        #         s_calling_cmd+=sf_tmp_cmd
        #     l_tmp_sh=gnrt_pipelines(s_head, s_libs, s_calling_cmd, sf_rep_sample_id, sf_rep_bam, sf_rep_x10_bam,
        #                    sf_sample_folder, rep_type)
        #     for tmp_sh in l_tmp_sh:
        #         l_sh.append(tmp_sh)
        # i_rep_type=get_flag_by_rep_type(rep_type)
        sf_config = sf_sample_folder+"config"

        #s_wfolder_rep=s_wfolder+rep_type
        # s_wfolder_rep = sf_sample_folder + rep_type
        if os.path.exists(sf_sample_folder)==False:
            cmd="mkdir {0}".format(sf_sample_folder)
            run_cmd(cmd)

        sf_rep_sample_id = sf_sample_folder + "/sample_id.txt"
        sf_rep_bam = sf_sample_folder + "/bam_list1.txt"
        if os.path.isfile(sf_rep_bam)==False:
            sf_rep_bam="null"
        sf_rep_x10_bam = sf_sample_folder + "/x10_bam_list1.txt"
        if os.path.exists(sf_rep_x10_bam)==False:
            sf_rep_x10_bam="null"
        
        s_head = ""
        if b_lsf==True:
            s_head = gnrt_script_head_lsf(spartition, ncores, stime, smemory, sid_tmp)
        elif b_slurm==True:
            s_head = gnrt_script_head(spartition, ncores, stime, smemory, sid_tmp, email_user)

        l_libs = load_par_config(sf_config)
        s_libs = gnrt_parameters(l_libs)
        ##
        iclip_c = args.nclip
        iclip_rp = args.cliprep
        idisc_c = args.ndisc
        icns_c = args.ncns # YW 2021/05/20 added
        # iflt_clip = args.nfilterclip
        # iflt_disc = args.nfilterdisc
        iflk_len = args.flklen
        # itei_len = args.teilen
        iflag = args.flag
        # YW 2021/05/26 added to do cns realignment in parallel
        c_realn_partition = args.c_realn_partition
        c_realn_time = args.c_realn_time
        c_realn_mem = args.c_realn_mem
        d_realn_partition = args.d_realn_partition
        d_realn_time = args.d_realn_time
        d_realn_mem = args.d_realn_mem
        check_interval = args.check_interval
        email_user = args.email_user
        # YW 2021/08/04 added to enable control bam file coordinate lifting
        sf_ref_bed = args.ref_bed
        error_margin = args.error_margin
        feat_extract_time = args.feat_extract_time
        feat_extract_partition = args.feat_extract_partition
        
        # YW 2021/05/20 removed b_mosaic, b_tumor, f_purity, i_rep_type, iflt_clip, iflt_disc, itei_len, added l_rep_type
        # YW 2021/08/04 added b_ctrl, sf_ref_bed, and error_margin
        s_calling_cmd = gnrt_calling_command(iclip_c, iclip_rp, idisc_c, icns_c, ncores, iflk_len,
                                             iflag, b_user_par, b_force, b_resume, l_rep_type, sf_pub_clip, sf_pub_clip_locus,
                                             c_realn_partition, c_realn_time, c_realn_mem, d_realn_partition, d_realn_time, d_realn_mem, check_interval, b_ctrl, sf_ref_bed, error_margin, sf_locus_file, feat_extract_time, feat_extract_partition, email_user)
        # if rep_type is REP_TYPE_SVA:
        #     s_calling_cmd = gnrt_calling_command(iclip_c, iclip_rp, idisc_c, icns_c, iflt_clip, iflt_disc, ncores, iflk_len,
        #                                          itei_len, iflag, b_mosaic, b_user_par, b_force, b_tumor, b_resume, f_purity, i_rep_type,sf_pub_clip, True, False)
        # # YW 2020/08/16 added the REP_TYPE_L1 statement
        # if rep_type is REP_TYPE_L1:
        #     s_calling_cmd = gnrt_calling_command(iclip_c, iclip_rp, idisc_c, icns_c, iflt_clip, iflt_disc, ncores, iflk_len,
        #                                          itei_len, iflag, b_mosaic, b_user_par, b_force, b_tumor, b_resume, f_purity, i_rep_type,sf_pub_clip, False, True)
        # if rep_type is REP_TYPE_HERV or rep_type is REP_TYPE_MSTA:
        #     s_calling_cmd = gnrt_calling_command_non_RNA(iclip_c, iclip_rp, idisc_c, icns_c, iflt_clip, iflt_disc, ncores, iflk_len,
        #                                          itei_len, iflag, b_mosaic, b_user_par, sf_pub_clip)
        # if rep_type is REP_TYPE_MIT:
        #     s_calling_cmd = gnrt_calling_command_MT(iclip_c, iclip_rp, idisc_c, icns_c, iflt_clip, iflt_disc, ncores, iflk_len,
        #                                          itei_len, iflag, b_mosaic, b_user_par, sf_pub_clip)
        #
        # if sf_case_control_bam_list != "null":
        #     sf_ctr_list=s_wfolder_rep + "/bam_list2.txt"
        #     sf_case_ctrl_list=s_wfolder_rep + "/bam_list3.txt"
        #     sf_tmp_cmd=gnrt_calling_command_somatic_case_control(iflt_clip, iflt_disc, ncores, sf_ctr_list,
        #                                                          iflk_len, sf_case_ctrl_list, iflag, i_rep_type)
        #     s_calling_cmd+=sf_tmp_cmd
        l_tmp_sh=gnrt_pipelines(s_head, s_libs, s_calling_cmd, sf_rep_sample_id, sf_rep_bam, sf_rep_x10_bam, sf_sample_folder)
        for tmp_sh in l_tmp_sh:
            l_sh.append(tmp_sh)
    with open(sf_submit_sh, "w") as fout_submit:
        fout_submit.write("#!/bin/bash\n\n")
        for s_sh in l_sh:
            if b_lsf==False:
                fout_submit.write("sbatch < "+s_sh+"\n")
            else:
                fout_submit.write("bsub < " + s_sh + "\n")
####

####Input:
# m_ids: sample id dictionary
# YW 2021/05/19 changed directory where sample_id.txt and bam_list*.txt are written to the general directory, not to each rep
def split_bam_list(m_ids, sf_bams, sf_x10_bams, s_wfolder, sf_case_control_bam_list="null"):
    #load in sf_bams
    m_bams = {}
    if sf_bams != "null":
        with open(sf_bams) as fin_bams:
            for line in fin_bams:
                fields = line.split()
                sid = fields[0]
                if sid not in m_ids:
                    continue
                s_bam = fields[1]
                if sid not in m_bams:
                    m_bams[sid]=[]
                m_bams[sid].append(s_bam)

    #load in sf_x10_bams
    m_bams_10X = {}
    if sf_bams_10X != "null":
        with open(sf_x10_bams) as fin_bams_10X:
            for line in fin_bams_10X:
                fields = line.split()
                sid = fields[0]
                if sid not in m_ids:
                    continue
                s_bam = fields[1]
                s_barcode_bam = fields[2]
                m_bams_10X[sid] = (s_bam, s_barcode_bam)
    
    #load in the control bams
    m_control={}
    m_case_control={}
    if sf_case_control_bam_list!="null":
        with open(sf_case_control_bam_list) as fin_cc_bam_list:
            for line in fin_cc_bam_list:
                fields=line.split()
                if len(fields)<3:
                    continue
                s_id=fields[0]
                m_control[s_id]=[]
                for sf_tmp in fields[2:]:
                    m_control[s_id].append(sf_tmp)
                m_case_control[s_id]=[]
                for sf_tmp in fields[1:]:
                    m_case_control[s_id].append(sf_tmp)

    for sid_tmp in m_ids:
        sf_sample_folder = s_wfolder + sid_tmp + "/"
        if os.path.exists(sf_sample_folder) == False:
            cmd = "mkdir {0}".format(sf_sample_folder)
            run_cmd(cmd)        
        # YW 2021/05/19 modified from before, to accomodate that we now have only one folder for each sample
        # need to generate two files for each sample
        sf_rep_sample_id = sf_sample_folder + "sample_id.txt"
        with open(sf_rep_sample_id, "w") as fout_rep_sample_id:
            fout_rep_sample_id.write(sid_tmp)
        if sid_tmp in m_bams:
            sf_rep_bam = sf_sample_folder + "bam_list1.txt"
            with open(sf_rep_bam, "w") as fout_rep_bams:
                for sf_tmp_bam in m_bams[sid_tmp]:
                    fout_rep_bams.write(sid_tmp+"\t"+sf_tmp_bam+"\n")
        if sid_tmp in m_control:
            sf_rep_bam2 = sf_sample_folder + "bam_list2.txt" #save the control file list
            with open(sf_rep_bam2, "w") as fout_rep_bams2:
                for sf_tmp_bam in m_control[sid_tmp]:
                    fout_rep_bams2.write(sf_tmp_bam+"\n")
        if sid_tmp in m_control:
            sf_rep_bam3 = sf_sample_folder + "bam_list3.txt" #save the control file list
            with open(sf_rep_bam3, "w") as fout_rep_bams3:
                fout_rep_bams3.write(sid_tmp)
                for sf_tmp_bam in m_case_control[sid_tmp]:
                    fout_rep_bams3.write("\t"+sf_tmp_bam)
                fout_rep_bams3.write("\n")
        if sid_tmp in m_bams_10X:
            sf_rep_x10_bam = sf_sample_folder + "x10_bam_list1.txt"
            with open(sf_rep_x10_bam, "w") as fout_rep_x10_bams:
                fout_rep_x10_bams.write(sid_tmp+"\t"+m_bams_10X[sid_tmp][0]+"\t"+m_bams_10X[sid_tmp][1]+"\n")


# YW 2021/08/04 added to enable control bam file coordinate lifting
def split_ctrl_bed_list(m_ids, sf_ctrl_bed, s_wfolder, m_ctrls):
    if sf_ctrl_bed != "null":
        with open(sf_ctrl_bed) as fin:
            for lines in fin:
                fields = lines.split()
                sid = fields[0]
                if sid not in m_ids:
                    continue
                s_ctrl = fields[1]
                if sid not in m_ctrls:
                    m_ctrls[sid] = s_ctrl
                else:
                    sys.exit(f"Error: There are duplicate ctrl bed input for {sid}! Exiting...")
        # for sid_tmp in m_ids:
        #     sf_sample_folder = s_wfolder + sid_tmp + "/"
        #     f_ctrl_bed = sf_sample_folder + "ctrl_bed.txt"
        #     with open(f_ctrl_bed, "w") as fout_ctrl_bed:
        #         fout_ctrl_bed.write(sid_tmp+"\t"+m_ctrls[sid_tmp]+"\n")


####run the pipelines
def run_pipeline(l_rep_type, sample_id, s_wfolder):
    for rep_type in l_rep_type:
        sf_sbatch_sh_rep = s_wfolder + sample_id + "/run_{0}.sh".format(rep_type)
        cmd="sh {0}".format(sf_sbatch_sh_rep)
        run_cmd(cmd)

def decompress(sf_in_tar, sf_out):
    cmd="tar -zxvf {0} -C {1}".format(sf_in_tar, sf_out)
    run_cmd(cmd)

def cp_compress_results(s_wfolder, l_rep_type, sample_id):
    # create a "results" folder
    sf_rslts = s_wfolder + "results/"
    if os.path.exists(sf_rslts)==False:
        cmd = "mkdir {0}".format(sf_rslts)
        run_cmd(cmd)

    for rep_type in l_rep_type:
        sf_rslts_rep_folder=sf_rslts+rep_type+"/"
        if os.path.exists(sf_rslts_rep_folder)==False:
            cmd = "mkdir {0}".format(sf_rslts_rep_folder)
            run_cmd(cmd)
        sf_samp_folder = sf_rslts_rep_folder + sample_id + "/"
        if os.path.exists(sf_samp_folder)==False:
            cmd="mkdir {0}".format(sf_samp_folder)
            run_cmd(cmd)

        sf_source_folder=s_wfolder+sample_id+"/"+rep_type+"/"
        sf_rslt1=sf_source_folder+"candidate_disc_filtered_cns.txt"
        cp_file(sf_rslt1, sf_samp_folder)
        sf_rslt2 = sf_source_folder + "candidate_list_from_clip.txt"
        cp_file(sf_rslt2, sf_samp_folder)
        sf_rslt2_1 = sf_source_folder + "candidate_list_from_clip.txt_tmp"
        cp_file(sf_rslt2_1, sf_samp_folder)
        sf_rslt3 = sf_source_folder + "candidate_list_from_disc.txt"
        cp_file(sf_rslt3, sf_samp_folder)

        # s_tmp1=sf_source_folder+"tmp/cns/candidate_sites_all_disc.fa"
        # cp_file(s_tmp1, sf_samp_folder)
        # s_tmp2 = sf_source_folder + "tmp/cns/candidate_sites_all_clip.fq"
        # cp_file(s_tmp2, sf_samp_folder)

        s_tmp1=sf_source_folder+"tmp/cns/temp_disc.sam"
        cp_file(s_tmp1, sf_samp_folder)
        s_tmp2 = sf_source_folder + "tmp/cns/temp_clip.sam"
        cp_file(s_tmp2, sf_samp_folder)
        s_tmp3 = sf_source_folder + "tmp/cns/all_with_polymerphic_flanks.fa"
        cp_file(s_tmp3, sf_samp_folder)
        s_tmp4 = sf_source_folder + "tmp/clip_reads_tmp0"
        cp_file(s_tmp4, sf_samp_folder)
        s_tmp5 = sf_source_folder + "tmp/discordant_reads_tmp0"
        cp_file(s_tmp5, sf_samp_folder)

    #compress the results folder to one file
    sf_compressed=sf_rslts+"results.tar.gz"
    cmd="tar -cvzf {0} -C {1} .".format(sf_compressed, sf_rslts)
    run_cmd(cmd)

def get_flag_by_rep_type(s_type):
    if s_type is REP_TYPE_L1:
        return 1
    elif s_type is REP_TYPE_ALU:
        return 2
    elif s_type is REP_TYPE_SVA:
        return 4
    elif s_type is REP_TYPE_HERV:
        return 8
    elif s_type is REP_TYPE_MIT:
        return 16
    elif s_type is REP_TYPE_MSTA:
        return 32

#
def prepare_case_control_bam(sf_ori_bam, sf_sprt_bam, sf_control_bam):
    with open(sf_ori_bam) as fin_ori, open(sf_sprt_bam, "w") as fout_sprt, open(sf_control_bam, "w") as fout_ctrl:
        for line in fin_ori:
            fields=line.split()
            if len(fields)<3:
                print("Not in right format. Less fields: {0}".format(line))
                continue
            s_id=fields[0]
            sf_case=fields[1]
            fout_sprt.write(s_id+"\t"+sf_case+"\n")

            fout_ctrl.write(s_id)
            for sf_tmp in fields[2:]:
                fout_ctrl.write("\t"+sf_tmp)
            fout_ctrl.write("\n")

####
def parse_arguments():
    parser = argparse.ArgumentParser("Generates xTEA sbatch scripts for hg38 to be run locally")
    parser.add_argument("-U", "--user",
                        action="store_true", dest="user", default=False,
                        help="Use user specific parameters instead of automatically calculated ones")
    parser.add_argument("--force",
                        action="store_true", dest="force", default=False,
                        help="Force to start from the very beginning")
    # parser.add_argument("--tumor",
    #                   action="store_true", dest="tumor", default=False,
    #                   help="Working on tumor samples")
    # parser.add_argument("--purity", dest="purity", type=float, default=0.45,  # by default tumor purity set to 45%
    #                   help="Tumor purity")
    parser.add_argument("--lsf",
                        action="store_true", dest="lsf", default=False,
                        help="Indicates submit to LSF system")
    parser.add_argument("--slurm",
                        action="store_true", dest="slurm", default=True,
                        help="Indicates submit to slurm system")
    parser.add_argument("--resume",
                        action="store_true", dest="resume", default=False,
                        help="Resume the running, which will skip the step if output file already exists!")
    parser.add_argument("-V", "--version",
                        action="store_true", dest="version", default=False,
                        help="Print xTea version")
    parser.add_argument("-i", "--id", dest="id",
                        help="sample id list file ", metavar="FILE")
    parser.add_argument("-a", "--par", dest="parameters",
                        help="parameter file ", metavar="FILE")
    parser.add_argument("-l", "--lib", dest="lib",
                        help="TE lib config file ", metavar="FILE")
    parser.add_argument("-b", "--bam", dest="bam",
                        help="Input bam file", metavar="FILE")
    parser.add_argument("-x", "--x10", dest="x10",
                        help="Input 10X bam file", metavar="FILE")
    parser.add_argument("-n", "--cores", dest="cores", type=int, default=8,
                        help="number of cores")
    parser.add_argument("-m", "--memory", dest="memory", type=int, default=20,
                        help="Memory limit in GB")
    parser.add_argument("-q", "--partition", dest="partition", type=str,
                        help="Which queue to run the job")
    parser.add_argument("-t", "--time", dest="time", type=str,
                        help="Time limit")

    parser.add_argument("-p", "--path", dest="wfolder", type=str, default="",
                        help="Working folder")
    parser.add_argument("-r", "--ref", dest="ref", type=str,
                        help="reference genome")
    parser.add_argument("-g", "--gene", dest="gene", type=str, default="gene.aff3",
                        help="Gene annotation file")
    parser.add_argument("--xtea", dest="xtea", type=str,
                        help="xTEA folder")

    parser.add_argument("-f", "--flag", dest="flag", type=int, default=5, # YW 2021/05/22 added default=3
                        # help="Flag indicates which step to run (1-clip, 2-disc, 4-barcode, 8-xfilter, 16-filter, 32-asm)")
                        help="Flag indicates which step to run (1-clip, 2-locus-specific-clip, 4-disc)")

    parser.add_argument("-y", "--reptype", dest="rep_type", type=int, default=7, # YW 2021/05/19 changed from 1
                        help="Type of repeats working on: 1-L1, 2-Alu, 4-SVA, 8-HERV, 16-Mitochondrial")

    parser.add_argument("--flklen", dest="flklen", type=int, default=3000,
                        help="flank region file")
    parser.add_argument("--nclip", dest="nclip", type=int, default=3,
                        help="cutoff of minimum # of clipped reads")
    parser.add_argument("--cr", dest="cliprep", type=int, default=1,
                        help="cutoff of minimum # of clipped reads whose mates map in repetitive regions")
    parser.add_argument("--nd", dest="ndisc", type=int, default=5,
                        help="cutoff of minimum # of discordant pair")
    # YW 2021/05/19 added the 2 arguments below
    parser.add_argument("--ncns", dest="ncns", type=int, default=1,
                        help="cutoff of minimum # of clipped reads + discordant pair of each sample aligned to repeat cns")
    parser.add_argument("--mat", "--final_matrix", dest="final_matrix",
                        help="The final feature matrix", metavar="FILE")
    # parser.add_argument("--nfclip", dest="nfilterclip", type=int, default=3,
    #                   help="cutoff of minimum # of clipped reads in filtering step")
    # parser.add_argument("--nfdisc", dest="nfilterdisc", type=int, default=5,
    #                   help="cutoff of minimum # of discordant pair of each sample in filtering step")
    # parser.add_argument("--teilen", dest="teilen", type=int, default=50,
    #                   help="minimum length of the insertion for future analysis")

    parser.add_argument("-o", "--output", dest="output", default="submit_calling_jobs_for_samples.sh",
                        help="The output file", metavar="FILE")
    parser.add_argument("--blacklist", dest="blacklist", default="null",
                        help="Reference panel database for filtering, or a blacklist region", metavar="FILE")
    
    # YW 2021/05/26 added to parallelize cns remapping
    parser.add_argument("--c_realn_partition", dest="c_realn_partition", type=str, default="short", # medium for highcov
                        help="slurm partition for running realignment of clipped reads to repeat cns")
    parser.add_argument("--c_realn_time", dest="c_realn_time", type=str, default="0-8:00", # 1-00:00 for highcov
                        help="runtime for running realignment of clipped reads to repeat cns")
    parser.add_argument("--c_realn_mem", dest="c_realn_mem", type=int, default=50,
                        help="run memory (in GB) for running realignment of clipped reads to repeat cns")
    parser.add_argument("--d_realn_partition", dest="d_realn_partition", type=str, default="short",
                        help="slurm partition for running realignment of disc reads to repeat cns")
    parser.add_argument("--d_realn_time", dest="d_realn_time", type=str, default="0-0:10",
                        help="runtime for running realignment of disc reads to repeat cns")
    parser.add_argument("--d_realn_mem", dest="d_realn_mem", type=int, default=1,
                        help="run memory (in GB) for running realignment of disc reads to repeat cns")
    parser.add_argument("--check_interval", dest="check_interval", type=int, default=60,
                        help="interval (in seconds) of checking whether the sbatch jobs of cns alignment of other repeat type have finished after finishing the main cns realignment job")
    # YW added 2023/03/29
    parser.add_argument("--email_user", dest="email_user", type=str, required=True,
                        help="specify the email address to notify SLURM job status")

    # YW 2021/08/04 added to enable control bam file coordinate lifting
    parser.add_argument("--ctrl", dest="ctrl",
                        action="store_true", default=False,
                        help="indicate running of aDNA ctrl bam input, skip feature extraction in -D step")
    parser.add_argument("--ref_bed", dest="ref_bed", default="/n/data1/bch/genetics/lee/elain/xTEA_benchmarking/alt_reference/gold_std_hg38_v4_int.temp",
                        help="The file with gold std set coordinate and insertion size, required only when --ctrl", metavar="FILE")
    parser.add_argument("--error_margin", dest="error_margin", type=str, default=15,
                        help="error margin to lift coordinate, required only when --ctrl")
    parser.add_argument("--ctrl_bed", dest="ctrl_bed", default="null",
                        help="File containing link to TEI coordinate from the control bam file, ancient only, required only when --ctrl==False", metavar="FILE")
    
    # YW 2021/10/12 added to enable extraction of features from predefined loci
    parser.add_argument("--locus_file", dest="locus_file", default="null",
                        help="The file with TEI loci to extract features, required only when -L/--locus_clip", metavar="FILE")
    
    # YW 2021/10/27 added to enable extraction of features without a matched reference control
    parser.add_argument("--train",
                        action="store_true", dest="train", default=False,
                        help="Turn on training mode, where there is a matched reference control")
    
    # # YW 2021/12/10 added to enable changing time limit for feature extraction by chunk
    parser.add_argument("--feat_extract_time", dest="feat_extract_time", default="0-01:30",
                        help="runtime for running feature extraction on 5M TEI sites")
    parser.add_argument("--feat_extract_partition", dest="feat_extract_partition", default="short",
                        help="partition for running feature extraction on 5M TEI sites")
    
    args = parser.parse_args()
    return args

####
if __name__ == '__main__':
    args = parse_arguments()
    b_version=args.version
    if b_version==True:
        print(("xTea %s for short and linked reads on hg38\n" % S_VERSION))
    else:
        sf_id = args.id
        sf_bams = args.bam ###input is a bam file
        sf_bams_10X = args.x10
        sf_sbatch_sh = args.output  # this is the shell for submitting the jobs
        s_wfolder = args.wfolder
        # b_mosaic = args.mosaic
        b_user_par = args.user
        b_force = args.force
        # b_case_control=args.case_control #whether case control mode
        # b_tumor = args.tumor
        # f_purity = args.purity
        b_lsf=args.lsf
        if b_lsf: # YW 2021/05/26 added because the cns realignment code is only written for SLURM
            raise NotImplementedError
        b_slurm=args.slurm ####
        b_resume = args.resume
        email_user = args.email_user
        
        # new args added by YW 2021/08/04 to enable control bam file coordinate lifting
        # YW 2021/10/28 added args.train if statement
        b_ctrl = args.ctrl
        sf_ctrl_bed = args.ctrl_bed
        if args.train:
            if not b_ctrl:
                if not sf_ctrl_bed:
                    sys.exit("TE calling results from control bam file is required!")
                if not os.path.isfile(sf_ctrl_bed):
                    sys.exit(f"{sf_ctrl_bed} doesn't exist. Exiting...")
            else:
                sf_ctrl_bed = "null"
            
        # new args added by YW 2021/10/12 to enable TEI feature collection from predefined loci
        sf_locus_file = args.locus_file
        
    ####

        if s_wfolder[-1]!="/":
            s_wfolder+="/"
        if os.path.exists(s_wfolder) == False:
            scmd = "mkdir {0}".format(s_wfolder)
            Popen(scmd, shell=True, stdout=PIPE).communicate()

        if os.path.isfile(sf_bams) == False:
            sf_bams = "null"
        if os.path.isfile(sf_bams_10X) == False:
            sf_bams_10X = "null"
    ####
        spartition = args.partition
        stime = args.time
        smemory = args.memory
        ncores = args.cores
        sf_folder_rep1 = args.lib  ##this is the lib folder path
        sf_ref1=args.ref ####reference genome
        sf_folder_xtea=args.xtea#

        sf_folder_rep=sf_folder_rep1
        sf_ref=sf_ref1
        # if args.decompress==True:
        #     decompress(sf_folder_rep1, s_wfolder)
        #     decompress(sf_ref1, s_wfolder)
        #     sf_folder_rep = s_wfolder+"rep_lib_annotation/" #trim tar.gz
        #     sf_ref=s_wfolder+"genome.fa"
        sf_gene=args.gene
        sf_black_list = args.blacklist

        i_rep_type=args.rep_type
        l_rep_type = []
        if i_rep_type & 1 != 0:
            l_rep_type.append(REP_TYPE_L1)
        if i_rep_type & 2 != 0:
            l_rep_type.append(REP_TYPE_ALU)
        if i_rep_type & 4 != 0:
            l_rep_type.append(REP_TYPE_SVA)
        if i_rep_type & 8 != 0:
            l_rep_type.append(REP_TYPE_HERV)
        if i_rep_type & 16 != 0:
            l_rep_type.append(REP_TYPE_MIT)
        if i_rep_type & 32 != 0:
            l_rep_type.append(REP_TYPE_MSTA)

        # if b_case_control==False:
        #     gnrt_running_shell(sf_id, sf_bams, sf_bams_10X, l_rep_type, b_mosaic, b_user_par, b_force, b_tumor, f_purity,
        #                        s_wfolder, sf_folder_rep, sf_ref, sf_gene, sf_black_list, sf_folder_xtea, spartition, stime,
        #                        smemory, ncores, sf_sbatch_sh, "null", b_lsf, b_slurm, b_resume)
        # else:
        #     #if in case-control mode, then need to separate to two steps
        #     #thus need to prepare the files separately
        #     sf_sprt_bams=sf_bams + CASE_LIST_SUFFIX
        #     sf_control_bams=sf_bams + CTRL_LIST_SUFFIX
        #     print(sf_control_bams)
        #     prepare_case_control_bam(sf_bams, sf_sprt_bams, sf_control_bams)
        #     #first, run the jobs seperately for all the case and control samples
        #     gnrt_running_shell(sf_id, sf_sprt_bams, sf_bams_10X, l_rep_type, b_mosaic, b_user_par, b_force, b_tumor,
        #                        f_purity, s_wfolder, sf_folder_rep, sf_ref, sf_gene, sf_black_list, sf_folder_xtea,
        #                        spartition, stime, smemory, ncores, sf_sbatch_sh, sf_bams, b_lsf, b_slurm, b_resume)
        gnrt_running_shell(sf_id, sf_bams, sf_bams_10X, l_rep_type, b_user_par, b_force,
                           s_wfolder, sf_folder_rep, sf_ref, sf_gene, sf_black_list, sf_folder_xtea, spartition, stime,
                           smemory, ncores, email_user, sf_sbatch_sh, "null", b_lsf, b_slurm, b_resume, b_ctrl, sf_ctrl_bed, sf_locus_file)

    #run_pipeline(l_rep_type, sample_id, s_wfolder)
    #cp_compress_results(s_wfolder, l_rep_type, sample_id)
####
####
####
