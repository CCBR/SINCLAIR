#!/usr/bin/env python

# GEX
# checks that all samples in INPUTDIR follow the format
# [Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz
# This script is based on the example at:
# https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv

# Testing
# bin/check_samplesheet.py assets/input_manifest.csv assets/contrast_manifest.csv /data/sevillas2/tmp/scRNA/project

import os
import sys
import errno
import argparse
import gzip
import re
from os import listdir
from os.path import isfile, join
import json
import pandas as pd
import numpy as np


def parse_args(args=None):
    Description = (
        "Samplesheets check - check both the samplesheet and contrast sheet contents."
    )
    Epilog = (
        "Example usage: python check_samplesheet.py <FILE_IN_S> <FILE_IN_C> <FILE_OUT>"
    )

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN_S", help="Input samplesheet file.")
    parser.add_argument("FILE_IN_C", help="Input contrast samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)


def check_files(fileid):
    # ensure that the file follows the structure
    ## [Sample Name]_S[Sample Number]_L00[Lane Number]_[Read Type]_001.fastq.gz
    # create sample list - will check all samples are the same name
    sample_list = list()
    sampleID = re.split("_S.", fileid)[0]
    sample_list.append(sampleID)

    # check S1_L00 is within file
    sID = re.search("S._L00", fileid)
    if sID == None:
        print_error("Input file must include S[sampleNumber]_L00 in name:", fileid)

    # spit after S1_L00
    postID = fileid.split("_L00")[1]

    # check lane, readype and extension
    laneID = postID.split("_")[0]
    try:
        float(laneID)
    except ValueError:
        print_error("Input file laneID must be numeric: %s" % fileid)

    accepted_read_type = ["R1", "R2", "I1", "I2"]
    read_type = postID.split("_")[1]
    if not read_type in accepted_read_type:
        print_error(
            "Input file read_type must be R1/R2/I1/I2 but %s was given" % read_type
        )

    ext = fileid.split("001.")[1]
    if ext != "fastq.gz":
        print_error("Input file extension must be .fastq.gz: %s" % ext)

    return sample_list


def remove_output_files(fname):
    if os.path.isfile(fname):
        os.remove(fname)


def check_samplesheet(file_in_s, file_in_c, file_out):
    """
    1) This function checks that the samplesheet follows the following structure:

    masterID,uniqueID,groupID,dataType,input_dir
    WB_Lysis_1,sample1,group1,gex,/data/CCBR_Pipeliner/Pipelines/TechDev_scRNASeq_Dev2023/test_dir/WB_Lysis_Granulocytes_3p_Introns_8kCells_fastqs/sample1
    WB_Lysis_1,sample2,group1,gex,/data/CCBR_Pipeliner/Pipelines/TechDev_scRNASeq_Dev2023/test_dir/WB_Lysis_Granulocytes_3p_Introns_8kCells_fastqs/sample2

    For an example see:
    https://github.com/CCBR/TechDev_scRNASeq_Dev2023/blob/dcd7f8e5b5ecabf48a3bb5b45ca73b8e0cde72c2/assets/input_manifest.csv

    2) This function checks that the contrast samplesheet follows the following structure:

    contrast1,contrast2
    WB_Lysis_1,WB_Lysis_2

    For an example see:
    https://github.com/CCBR/TechDev_scRNASeq_Dev2023/blob/dcd7f8e5b5ecabf48a3bb5b45ca73b8e0cde72c2/assets/contrast_manifest.csv

    Example usage:
    python bin/check_samplesheet.py assets/input_manifest.csv assets/contrast_manifest.csv /data/sevillas2/tmp/project; \
    python bin/check_samplesheet.py assets/input_manifest_notunique.csv assets/contrast_manifest.csv /data/sevillas2/tmp/project
    """

    ###################################################################################################
    # ALL: Check sample manifest
    ###################################################################################################
    mani_mapping_dict = {}
    group_mapping_dict = {}
    with open(file_in_s, "r") as fin:
        ## Check header
        MIN_COLS = 2
        HEADER = ["masterID", "uniqueID", "groupID", "dataType", "input_dir"]
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        if header[: len(HEADER)] != HEADER:
            print(
                "ERROR: Please check samplesheet header -> {} != {}".format(
                    ",".join(header), ",".join(HEADER)
                )
            )
            sys.exit(1)

        ## Check sample entries
        for line in fin.readlines():
            lspl = [x.strip().strip('"') for x in line.strip().split(",")]

            # If it's a blank line, next
            if len(line.strip()) == 0:
                print("Skipping blank line")
                continue

            # Check valid number of columns per row
            if len(lspl) < len(HEADER):
                print_error(
                    "Invalid number of columns (minimum = {})!".format(len(HEADER)),
                    "Line",
                    line,
                )
            num_cols = len([x for x in lspl if x])
            if num_cols < MIN_COLS:
                print_error(
                    "Invalid number of populated columns (minimum = {})!".format(
                        MIN_COLS
                    ),
                    "Line",
                    line,
                )

            # split lines
            MASTERID, UNIQUEID, GROUPID, DATATYPE, INPUTDIR = lspl[: len(HEADER)]

            ## Check sample name entries
            UNIQUEID = UNIQUEID.replace(" ", "_")
            if not UNIQUEID:
                print_error("Sample entry has not been specified!", "Line", line)

            ## Check input dir exists, run tests for files
            if os.path.exists(INPUTDIR):
                # check all files in input
                onlyfiles = [f for f in listdir(INPUTDIR) if isfile(join(INPUTDIR, f))]

                ###################################################################################################
                # GEX: check manifest
                ###################################################################################################
                if DATATYPE == "gex":
                    # for every file, check file exists and is properly formatted
                    for fileid in onlyfiles:
                        sample_list = check_files(fileid)

                    # ensure all samples names are the same
                    if len(set(sample_list)) != 1:
                        print_error(
                            "All input files within given dir must have the same sample ID %s"
                            % sample_list
                        )

                ###################################################################################################
                # ATAC: check manifest
                ###################################################################################################
                if DATATYPE == "atac":
                    # for every tile
                    sample_list = list()
                    for fileid in onlyfiles:
                        # ensure that the file follows the structure
                        ## [Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz
                        # create sample list - will check all samples are the same name
                        sampleID = fileid.split("_S1")[0]
                        sample_list.append(sampleID)

                    # ensure all samples names are the same
                    if len(set(sample_list)) != 1:
                        print_error(
                            "All input files within given dir must have the same sample ID %s"
                            % sample_list
                        )

            ###################################################################################################
            # ALL: Create DICTs for
            ## 1) MANIFEST INPUT
            ## 2) GROUP info
            ###################################################################################################
            ## Create manifest mapping dictionary
            ### uniqueid: data_tyoe: [INPUTDIR]
            if UNIQUEID not in mani_mapping_dict:
                mani_mapping_dict[UNIQUEID] = {}
            else:
                print_error(
                    "Samplesheet contains duplicate sample names!", "Line", line
                )
            mani_mapping_dict[UNIQUEID][DATATYPE] = {}
            mani_mapping_dict[UNIQUEID][DATATYPE] = [INPUTDIR]

            ## Create sample mapping dictionary
            ### groupid: data_tyoe: [uniqueid]
            if GROUPID not in group_mapping_dict:
                group_mapping_dict[GROUPID] = {}
            if DATATYPE not in group_mapping_dict[GROUPID]:
                group_mapping_dict[GROUPID][DATATYPE] = {}
            if MASTERID not in group_mapping_dict[GROUPID][DATATYPE]:
                group_mapping_dict[GROUPID][DATATYPE][MASTERID] = []
            try:
                group_mapping_dict[GROUPID][DATATYPE][MASTERID].append(UNIQUEID)
            except KeyError:
                group_mapping_dict[GROUPID][DATATYPE][MASTERID] = [UNIQUEID]

    ###################################################################################################
    # ALL: Create DICTs for
    ## 1) Group level sample INPUT
    ###################################################################################################
    ## for each of the groupIDS
    ## find the GEX data-types
    ## for each masterID, generate list of samples matching
    ## NOTE: samples can only be added if there is a GEX sample associated
    ## IE if there is an ATAC sample with masterID1, but no GEX sample, this sample will not be included
    group_sample_mapping_dict = {}
    for GROUPID in sorted(group_mapping_dict.keys()):
        for DATATYPE in sorted(group_mapping_dict[GROUPID].keys()):
            if DATATYPE == "gex":
                for MASTERID in sorted(group_mapping_dict[GROUPID][DATATYPE].keys()):
                    for UNIQUEID in group_mapping_dict[GROUPID][DATATYPE][MASTERID]:
                        if GROUPID not in group_sample_mapping_dict:
                            group_sample_mapping_dict[GROUPID] = []
                        try:
                            group_sample_mapping_dict[GROUPID].append(UNIQUEID)
                        except KeyError:
                            group_sample_mapping_dict[GROUPID] = [UNIQUEID]

    # ###################################################################################################
    # # ALL: Check manifest
    # ###################################################################################################
    contrast_mapping_dict = {}
    with open(file_in_c, "r") as fin:
        ## Check header
        MIN_COLS = 2
        header = [x.strip('"') for x in fin.readline().strip().split(",")]

        ## Check sample entries
        for line in fin:
            lspl = [x.strip().strip('"') for x in line.strip().split(",")]

            # Check valid number of columns per row
            num_cols = len([x for x in lspl if x])
            if num_cols < MIN_COLS:
                print_error(
                    "Invalid number of populated columns (minimum = {})!".format(
                        MIN_COLS
                    ),
                    "Line",
                    line,
                )

            for cid in lspl:
                if cid not in group_mapping_dict:
                    print_error(
                        "CONTRAST entry is not listed in the sample manifest. Check names!",
                        "Line",
                        line,
                    )

        ###################################################################################################
        # ALL: Create Contrast DF
        ###################################################################################################
        contrast_df = pd.read_csv(file_in_c)
        cols_to_cat = contrast_df.columns.values.tolist()
        contrast_df["key"] = contrast_df.astype(str).agg("-".join, axis=1)
        contrast_df["key"] = contrast_df["key"].str.replace("-nan", "")

    ###################################################################################################
    # ALL: Create OUTPUT
    ###################################################################################################
    # cleanup outputs if they exist
    remove_output_files(file_out + "_gex_samplesheet.csv")
    remove_output_files(file_out + "_contrast_samplesheet.csv")
    remove_output_files(file_out + "_atac_samplesheet.csv")
    remove_output_files(file_out + "*_groups_samplesheet.csv")

    # Write validated samplesheet with appropriate columns for cellranger
    for sample in sorted(mani_mapping_dict.keys()):
        for dt in sorted(mani_mapping_dict[sample].keys()):
            # GEX samplesheet
            if dt == "gex":
                fname = file_out + "_gex_samplesheet.csv"
                if not os.path.isfile(fname):
                    with open(fname, "w") as fout:
                        fout.write(",".join(["sample", "gex_input_dir"]) + "\n")
                        fout.close()
                with open(fname, "a+") as fout:
                    for idir in mani_mapping_dict[sample][dt]:
                        fout.write(sample + "," + idir + "\n")
                        fout.close()
            # ATAC samplesheet
            if dt == "atac":
                fname = file_out + "_atac_samplesheet.csv"
                if not os.path.isfile(fname):
                    with open(fname, "w") as fout:
                        fout.write(",".join(["sample", "atac_input_dir"]) + "\n")
                        fout.close()
                with open(fname, "a+") as fout:
                    for idir in mani_mapping_dict[sample][dt]:
                        fout.write(sample + "," + idir + "\n")
                        fout.close()

    # Write validated contrast samplesheet
    fname = file_out + "_contrast_samplesheet.csv"
    if not os.path.isfile(fname):
        contrast_df.to_csv(fname, sep=",", index=False)

    # Write validated group samplesheets
    fname = file_out + "_groups_samplesheet.csv"
    with open(fname, "w") as fout:
        fout.write(",".join(["keyid", "sampleid"]) + "\n")
        fout.close()
    for keyid in contrast_df["key"]:
        with open(fname, "a+") as fout:
            gid_list = keyid.split("-")
            for gid in gid_list:
                for sid in group_sample_mapping_dict[gid]:
                    fout.write(keyid + "," + sid + "\n")
    fout.close()


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN_S, args.FILE_IN_C, args.FILE_OUT)


if __name__ == "__main__":
    sys.exit(main())
