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
S_VERSION="v0.1"

####
PUB_CLIP="pub_clip"
PUB_CLIP_LOCUS="pub_clip_locus"
REP_TYPE_L1="L1"
REP_TYPE_ALU="Alu"
REP_TYPE_SVA="SVA"

# YW 2020/07/24 customized memory usage by rep_type
# YW 2021/05/20 removed rep_type
def gnrt_script_head(spartition, ncores, stime, smemory, s_id, email_user, s_wfolder):
    s_head = "#!/bin/bash\n\n"
    s_head += f"#SBATCH -c {ncores}\n"
    s_head += f"#SBATCH -t {stime}\n"
    s_head += f"#SBATCH --mem={smemory}G\n"
    s_head += f"#SBATCH -p {spartition}\n"
    s_head += f"#SBATCH -J {s_id}\n" # YW added 2021/05/24
    s_head += f"#SBATCH -o {s_wfolder}{s_id}/%j.out\n"
    s_head += "#SBATCH --mail-type=END,FAIL\n"
    s_head += f"#SBATCH --mail-user={email_user}\n"
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

def run_cmd(cmd):
    print(cmd)
    Popen(cmd, shell=True, stdout=PIPE).communicate()
    
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

    s_resume = ""
    if b_resume == True:
        s_resume = "--resume"
    
    setup_step = "cd ${PREFIX}\n"
    setup_step += "ln -nsf ${XTEA_PATH} ${TMP}source_scripts\n" # YW 2021/05/26 added to import modules from working directory
    sclip_step = f"python3 ${{XTEA_PATH}}\"x_TEA_main.py\" -C -i ${{BAM_LIST}} --lc {iclip_c} --rc {iclip_c} --cr {iclip_rp} " \
                 f"-o ${{PREFIX}}\"candidate_list_from_clip.txt\"  -n {ncores} --cp {s_cfolder} {s_user} {s_clean} {s_resume} " \
                 f"--ref ${{REF}} -p ${{TMP}} --c_realn_partition {c_realn_partition} --c_realn_time {c_realn_time} --c_realn_mem {c_realn_mem} --check_interval {check_interval} --email_user {email_user} --sample_id ${{SAMPLE_ID}} "
    if REP_TYPE_ALU in l_rep_type:
        sclip_step += "--Alu-annotation ${ALU_ANNOTATION} --Alu-reference ${ALU_COPY_WITH_FLANK} --Alu-cns ${ALU_CNS} "
    if REP_TYPE_L1 in l_rep_type:
        sclip_step += "--L1-annotation ${L1_ANNOTATION} --L1-reference ${L1_COPY_WITH_FLANK} --L1-cns ${L1_CNS} "
    if REP_TYPE_SVA in l_rep_type:
        sclip_step += "--SVA-annotation ${SVA_ANNOTATION} --SVA-reference ${SVA_COPY_WITH_FLANK} --SVA-cns ${SVA_CNS}\n"

    sdisc_step = f"python3 ${{XTEA_PATH}}\"x_TEA_main.py\" -D -i ${{PREFIX}}\"candidate_list_from_clip.txt\" --nd {idisc_c} " \
                 f"--ncns {icns_c} --ref ${{REF}} -b ${{BAM_LIST}} -p ${{TMP}} -o ${{PREFIX}}\"candidate_list_from_disc.txt\" -n {ncores} -m ${{PREFIX}}\"feature_matrix.txt\" {s_user} {s_resume} " \
                 f"--d_realn_partition {d_realn_partition} --d_realn_time {d_realn_time} --d_realn_mem {d_realn_mem} --check_interval {check_interval} --error_margin {error_margin} --email_user {email_user} --sample_id ${{SAMPLE_ID}} " # YW 2021/08/04 added error_margin
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
                 f"--ref ${{REF}} -p ${{TMP}} --c_realn_partition {c_realn_partition} --c_realn_time {c_realn_time} --c_realn_mem {c_realn_mem} --check_interval {check_interval} --email_user {email_user} --sample_id ${{SAMPLE_ID}} "
    if REP_TYPE_ALU in l_rep_type:
        sclip_locus_step += "--Alu-annotation ${ALU_ANNOTATION} --Alu-reference ${ALU_COPY_WITH_FLANK} --Alu-cns ${ALU_CNS} "
    if REP_TYPE_L1 in l_rep_type:
        sclip_locus_step += "--L1-annotation ${L1_ANNOTATION} --L1-reference ${L1_COPY_WITH_FLANK} --L1-cns ${L1_CNS} "
    if REP_TYPE_SVA in l_rep_type:
        sclip_locus_step += "--SVA-annotation ${SVA_ANNOTATION} --SVA-reference ${SVA_COPY_WITH_FLANK} --SVA-cns ${SVA_CNS}\n"
    
    sdisc_locus_step = f"python3 ${{XTEA_PATH}}\"x_TEA_main.py\" -D -i ${{PREFIX}}\"candidate_list_from_clip_locus.txt\" --nd {idisc_c} " \
                       f"--ncns {icns_c} --ref ${{REF}} -b ${{BAM_LIST}} -p ${{TMP}} -o ${{PREFIX}}\"candidate_list_from_disc_locus.txt\" -n {ncores} -m ${{PREFIX}}\"feature_matrix_locus.txt\" {s_user} {s_resume} " \
                       f"--d_realn_partition {d_realn_partition} --d_realn_time {d_realn_time} --d_realn_mem {d_realn_mem} --check_interval {check_interval} --error_margin {error_margin} --email_user {email_user} --sample_id ${{SAMPLE_ID}} "
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
    return s_cmd


####gnrt the whole pipeline
# YW 2021/05/23 removed rep_type, use the sf_root_folder as it is, not sf_working_folder, commented out code to generate cns, tsdct, igv folders
def gnrt_pipelines(s_head, s_libs, s_calling_cmd, sf_id, sf_bams, sf_root_folder):
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
            # if sid in m_bams_10X:
            #     fout_bam_list.write(m_bams_10X[sid] + "\t"+ X10 + "\n")
        ####gnrt the pipeline file
        sf_out_sh = sf_folder + "run_xTEA_pipeline.sh"
        with open(sf_out_sh, "w") as fout_sh:  ###gnrt the pipeline file
            fout_sh.write(s_head)
            s_prefix = f"PREFIX={sf_folder}\n"
            fout_sh.write(s_prefix)
            s_sid = f"SAMPLE_ID={sid}\n"
            fout_sh.write(s_sid)
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
def write_to_config(l_sf_anno, sf_ref, l_sf_copy_with_flank, sf_flank, l_sf_cns,
                    sf_xtea, s_bl, s_tmp, s_tmp_clip, sf_ctrl_bed, sf_config):
    with open(sf_config, "w") as fout:
        for sf_anno in l_sf_anno:
            fout.write(sf_anno)
        fout.write(sf_ref)
        for sf_copy_with_flank in l_sf_copy_with_flank:
            fout.write(sf_copy_with_flank)
        fout.write(sf_flank)
        for sf_cns in l_sf_cns:
            fout.write(sf_cns)
        fout.write(sf_xtea)
        fout.write(s_bl)
        fout.write(s_tmp)
        fout.write(s_tmp_clip)
        fout.write(sf_ctrl_bed)
#
#generate library configuration files
# YW 2021/05/19 enabled adding multiple library configuration files to the same command (separated by repeat type/MT)
# YW 2021/08/04 added sf_ctrl_bed to enable control bam file coordinate lifting
def gnrt_lib_config(l_rep_type, sf_folder_rep, sf_ref,  sf_folder_xtea, sf_config_prefix, sf_ctrl_bed):
    if sf_folder_rep[-1] != "/":
        sf_folder_rep += "/"
    if sf_folder_xtea[-1] != "/":
        sf_folder_xtea += "/"
    if sf_config_prefix[-1] != "/":
        sf_config_prefix += "/"

    s_bl = "BAM_LIST ${PREFIX}\"bam_list.txt\"\n"
    s_tmp = "TMP ${PREFIX}\"tmp/\"\n"
    s_tmp_clip = "TMP_CLIP ${PREFIX}\"tmp/clip/\"\n"
    sf_ref = "REF " + sf_ref + "\n"
    sf_xtea = "XTEA_PATH " + sf_folder_xtea + "\n"
    sf_flank = "SF_FLANK null\n"
    sf_ctrl_bed = "CTRL_BED " + sf_ctrl_bed + "\n"
    l_sf_anno = []
    # l_sf_anno1 = []
    l_sf_copy_with_flank = []
    l_sf_cns = []
    for s_rep_type in l_rep_type:
        if s_rep_type=="Alu":# for Alu
            l_sf_anno.append("ALU_ANNOTATION " + sf_folder_rep + "Alu/hg38/hg38_Alu.out\n")
            l_sf_copy_with_flank.append("ALU_COPY_WITH_FLANK " + sf_folder_rep + "Alu/hg38/hg38_AluJabc_copies_with_flank.fa\n")
            l_sf_cns.append("ALU_CNS " + sf_folder_rep + "consensus/ALU.fa\n")

        elif s_rep_type=="L1":# for L1
            l_sf_anno.append("L1_ANNOTATION " + sf_folder_rep + "L1/hg38/hg38_L1_larger_500_with_all_L1HS.out\n")
            l_sf_copy_with_flank.append("L1_COPY_WITH_FLANK " + sf_folder_rep + "L1/hg38/hg38_L1_copies_larger_5K_with_flank.fa\n")
            l_sf_cns.append("L1_CNS " + sf_folder_rep + "consensus/LINE1.fa\n")

        elif s_rep_type=="SVA":# for SVA
            l_sf_anno.append("SVA_ANNOTATION " + sf_folder_rep + "SVA/hg38/hg38_SVA.out\n")
            l_sf_copy_with_flank.append("SVA_COPY_WITH_FLANK " + sf_folder_rep + "SVA/hg38/hg38_SVA_copies_with_flank.fa\n")
            l_sf_cns.append("SVA_CNS " + sf_folder_rep + "consensus/SVA.fa\n")

    # YW 2021/05/20
    sf_config = sf_config_prefix + "config"
    write_to_config(l_sf_anno, sf_ref, l_sf_copy_with_flank, sf_flank,
                    l_sf_cns, sf_xtea, s_bl, s_tmp, s_tmp_clip, sf_ctrl_bed, sf_config) # YW 2021/08/04 added sf_ctrl_bed


####gnrt the running shell
# YW 2021/05/20 modified to generate calling command using one config file, not by different rep_type, removed b_mosaic, b_tumor, f_purity
# YW 2021/08/04 added b_ctrl=False, sf_ctrl_bed="null" to enable control bam file coordinate lifting
def gnrt_running_shell(sf_ids, sf_bams, l_rep_type, b_user_par, b_force, s_wfolder,
                       sf_folder_rep, sf_ref, sf_folder_xtea, spartition, stime, smemory,
                       ncores, email_user, sf_submit_sh, sf_case_control_bam_list="null", b_slurm=True, b_resume=False,
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
    split_bam_list(m_id, sf_bams, s_wfolder)
    
    # YW 2021/08/04 added to enable control bam file coordinate lifting, generate the ref_bed file
    m_ctrls = split_ctrl_bed_list(m_id, sf_ctrl_bed, s_wfolder)
    
    l_sh=[]
    for sid_tmp in m_id:
        sf_sample_folder=s_wfolder + sid_tmp + "/"
        sf_pub_clip = sf_sample_folder + PUB_CLIP + "/"
        sf_pub_clip_locus = sf_sample_folder + PUB_CLIP_LOCUS + "/"
        # YW 2021/08/04 added m_ctrls to enable control bam file coordinate lifting
        if sid_tmp in m_ctrls:
            gnrt_lib_config(l_rep_type, sf_folder_rep, sf_ref, sf_folder_xtea, sf_sample_folder, m_ctrls[sid_tmp])
        else:
            gnrt_lib_config(l_rep_type, sf_folder_rep, sf_ref, sf_folder_xtea, sf_sample_folder, "null")

        sf_config = sf_sample_folder+"config"

        if os.path.exists(sf_sample_folder)==False:
            cmd="mkdir {0}".format(sf_sample_folder)
            run_cmd(cmd)

        sf_rep_sample_id = sf_sample_folder + "/sample_id.txt"
        sf_rep_bam = sf_sample_folder + "/bam_list1.txt"
        if os.path.isfile(sf_rep_bam)==False:
            sf_rep_bam="null"
        s_head = gnrt_script_head(spartition, ncores, stime, smemory, sid_tmp, email_user, s_wfolder)

        l_libs = load_par_config(sf_config)
        s_libs = gnrt_parameters(l_libs)
        ##
        iclip_c = args.nclip
        iclip_rp = args.cliprep
        idisc_c = args.ndisc
        icns_c = args.ncns # YW 2021/05/20 added
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
        l_tmp_sh=gnrt_pipelines(s_head, s_libs, s_calling_cmd, sf_rep_sample_id, sf_rep_bam, sf_sample_folder)
        for tmp_sh in l_tmp_sh:
            l_sh.append(tmp_sh)
    with open(sf_submit_sh, "w") as fout_submit:
        fout_submit.write("#!/bin/bash\n\n")
        for s_sh in l_sh:
            fout_submit.write("sbatch < "+s_sh+"\n")
####

####Input:
# m_ids: sample id dictionary
# YW 2021/05/19 changed directory where sample_id.txt and bam_list*.txt are written to the general directory, not to each rep
def split_bam_list(m_ids, sf_bams, s_wfolder):
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


# YW 2021/08/04 added to enable control bam file coordinate lifting
def split_ctrl_bed_list(m_ids, sf_ctrl_bed, s_wfolder):
    m_ctrls = dict()
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
    return m_ctrls


####
def parse_arguments():
    parser = argparse.ArgumentParser("Generates xTEA sbatch scripts for hg38 to be run locally")
    parser.add_argument("-U", "--user",
                        action="store_true", dest="user", default=False,
                        help="Use user specific parameters instead of automatically calculated ones")
    parser.add_argument("--force",
                        action="store_true", dest="force", default=False,
                        help="Force to start from the very beginning")
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
    parser.add_argument("--xtea", dest="xtea", type=str,
                        help="xTEA folder")

    parser.add_argument("-f", "--flag", dest="flag", type=int, default=5, # YW 2021/05/22 added default=3
                        # help="Flag indicates which step to run (1-clip, 2-disc, 4-barcode, 8-xfilter, 16-filter, 32-asm)")
                        help="Flag indicates which step to run (1-clip, 2-locus-specific-clip, 4-disc)")

    parser.add_argument("-y", "--reptype", dest="rep_type", type=int, default=7, # YW 2021/05/19 changed from 1
                        help="Type of repeats working on: 1-L1, 2-Alu, 4-SVA")

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

    parser.add_argument("-o", "--output", dest="output", default="submit_calling_jobs_for_samples.sh",
                        help="The output file", metavar="FILE")
    # parser.add_argument("--blacklist", dest="blacklist", default="null",
    #                     help="Reference panel database for filtering, or a blacklist region", metavar="FILE")
    
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
        sf_sbatch_sh = args.output  # this is the shell for submitting the jobs
        s_wfolder = args.wfolder
        b_user_par = args.user
        b_force = args.force
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


        i_rep_type=args.rep_type
        l_rep_type = []
        if i_rep_type & 1 != 0:
            l_rep_type.append(REP_TYPE_L1)
        if i_rep_type & 2 != 0:
            l_rep_type.append(REP_TYPE_ALU)
        if i_rep_type & 4 != 0:
            l_rep_type.append(REP_TYPE_SVA)

        gnrt_running_shell(sf_id, sf_bams, l_rep_type, b_user_par, b_force,
                           s_wfolder, sf_folder_rep, sf_ref, sf_folder_xtea, spartition, stime,
                           smemory, ncores, email_user, sf_sbatch_sh, "null", b_slurm, b_resume, b_ctrl, sf_ctrl_bed, sf_locus_file)

####
####
####
