#!/usr/bin/env python

"""
Yleaf detection of Y-Haplogroups in Human DNA v3.1

Developed at Erasmus Medical Center Department of Genetic Identification

License: GNU General Public License v3 or later
A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

Autor: Diego Montiel Gonzalez
Extensively modified by: Bram van Wersch
"""

import argparse
import datetime
import io
import logging
import multiprocessing
import os
import re
import shutil
import subprocess
import sys
import time
from argparse import ArgumentParser
from collections import Counter, defaultdict
from functools import partial
from pathlib import Path
from typing import TextIO, Any

import numpy as np
import pandas as pd

from yleaf import __version__, download_reference, yleaf_constants
from yleaf.configuration import Configuration
from yleaf.tree import Tree
from yleaf.exceptions import YleafError, ExternalCommandError, ConfigurationError
from yleaf.utils import run_command

# pd.options.mode.chained_assignment = None  # Removed to fix issues properly

ACCEPTED_REF_BASES: set[str] = {"A", "C", "G", "T"}

START_RE = re.compile(r"\^.")
INDEL_RE = re.compile(r"([+-])(\d+)")

# path constants
PREDICTION_OUT_FILE_NAME: str = "hg_prediction.hg"
HAPLOGROUP_IMAGE_FILE_NAME: str = "hg_tree_image"

LOG: logging.Logger = logging.getLogger("yleaf_logger")


class MyFormatter(logging.Formatter):
    """
    Copied from MultiGeneBlast (my code)
    """

    def __init__(self, fmt, starttime=time.time()):
        logging.Formatter.__init__(self, fmt)
        self._start_time = starttime

    def format(self, record):
        """
        Overwrite of the format function that prints the passed time and adds
        current time to the existing format
        :See: logging.Formatter.format()
        """
        record.passedTime = f"{time.time() - self._start_time:.3f}"
        record.currentTime = datetime.datetime.now().time()
        return super(MyFormatter, self).format(record)


class YleafPipeline:
    def __init__(self, config: Configuration):
        self.config = config
        self.cached_position_data: set[str] | None = None
        self.cached_snp_database: dict[str, list[dict[str, str]]] | None = None
        self.cached_reference_file: str | None = None

    def run_vcf(
        self,
        path_markerfile: Path,
        base_out_folder: Path,
        args: argparse.Namespace,
        sample_vcf_file: Path,
    ):
        LOG.debug("Starting with extracting haplogroups...")
        markerfile = pd.read_csv(
            path_markerfile,
            header=None,
            sep="\t",
            names=["chr", "marker_name", "haplogroup", "pos", "mutation", "anc", "der"],
            dtype={
                "chr": str,
                "marker_name": str,
                "haplogroup": str,
                "pos": int,
                "mutation": str,
                "anc": str,
                "der": str,
            },
        ).drop_duplicates(subset="pos", keep="first", inplace=False)

        sample_vcf_folder = base_out_folder / (sample_vcf_file.name.replace(".vcf.gz", ""))
        self.safe_create_dir(sample_vcf_folder, args.force)

        cmd = ["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n", str(sample_vcf_file)]

        # Stream output to pandas
        # To avoid deadlock when reading from stdout while stderr fills up, we redirect stderr to a temp file.
        # This is robust for large outputs and large error logs.
        stderr_file = sample_vcf_folder / "bcftools.stderr"
        try:
            with open(stderr_file, "w") as err:
                with subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=err, text=True) as process:
                    try:
                        pileupfile = pd.read_csv(process.stdout, dtype=str, header=None, sep="\t")
                    except pd.errors.EmptyDataError:
                        pileupfile = pd.DataFrame()

                    process.wait()
                    if process.returncode != 0:
                        # Read the error file content for the exception
                        with open(stderr_file, "r") as err_read:
                            stderr_content = err_read.read()
                        raise ExternalCommandError(" ".join(cmd), process.returncode, stderr_content)
        except Exception as e:
            LOG.error(f"Failed to run bcftools query: {e}")
            raise ExternalCommandError(" ".join(cmd), -1, str(e))
        finally:
            if stderr_file.exists():
                os.remove(stderr_file)

        if pileupfile.empty:
             LOG.warning(f"No data returned from bcftools for {sample_vcf_file}")
             # Create placeholder output files to ensure downstream processing handles this sample as "NA"
             # instead of skipping it entirely.
             outputfile = sample_vcf_folder / (sample_vcf_file.name.replace(".vcf.gz", ".out"))
             fmf_output = sample_vcf_folder / (sample_vcf_file.name.replace(".vcf.gz", ".fmf"))

             # Create empty .out file with header
             with open(outputfile, "w") as f:
                 f.write(
                    "\t".join(
                        [
                            "chr",
                            "pos",
                            "marker_name",
                            "haplogroup",
                            "mutation",
                            "anc",
                            "der",
                            "reads",
                            "called_perc",
                            "called_base",
                            "state",
                            "depth\n",
                        ]
                    )
                )

             # Create empty .fmf file with header
             fmf_header = [
                "chr",
                "pos",
                "marker_name",
                "haplogroup",
                "mutation",
                "anc",
                "der",
                "reads",
                "called_perc",
                "called_base",
                "state",
                "Description",
            ]
             pd.DataFrame(columns=fmf_header).to_csv(fmf_output, sep="\t", index=False)

             # Write info file
             general_info_list = ["Total of mapped reads: VCF", "Total of unmapped reads: VCF"]
             general_info_list += ["Valid markers: 0"]
             general_info_list.append("Markers with zero reads: 0")
             general_info_list.append("Markers below the read threshold {0}: 0") # Thresholds not readily available here without args access or hardcoding default logic
             general_info_list.append("Markers below the base majority threshold {0}: 0")
             general_info_list.append("Markers with discordant genotype: 0")
             general_info_list.append("Markers without haplogroup information: 0")
             general_info_list.append("Markers with haplogroup information: 0")

             self.write_info_file(sample_vcf_folder, general_info_list)
             return

        pileupfile.columns = ["chr", "pos", "refbase", "altbase", "reads"]
        pileupfile["pos"] = pileupfile["pos"].astype(int)

        pileupfile["altbase"] = pileupfile["altbase"].str.split(",")
        pileupfile["reads"] = pileupfile["reads"].str.split(",")

        pileupfile["ref_reads"] = pileupfile["reads"].apply(lambda x: x[0])
        pileupfile["alt_reads"] = pileupfile["reads"].apply(lambda x: x[1:])

        pileupfile["alt_reads_dict"] = pileupfile.apply(
            lambda row: dict(zip(row["altbase"], row["alt_reads"])), axis=1
        )
        pileupfile["alt_reads_dict"] = pileupfile["alt_reads_dict"].apply(
            lambda x: {k: int(v) for k, v in x.items()}
        )
        pileupfile["highest_alt_reads"] = pileupfile["alt_reads_dict"].apply(
            lambda x: max(x.values()) if len(x) > 0 else 0
        )
        pileupfile["highest_alt_reads_base"] = pileupfile["alt_reads_dict"].apply(
            lambda x: max(x, key=x.get) if len(x) > 0 else "NA"
        )
        pileupfile["total_reads"] = pileupfile.apply(
            lambda row: int(row["ref_reads"]) + row["highest_alt_reads"], axis=1
        )
        pileupfile["called_ref_perc"] = pileupfile.apply(
            lambda row: round((int(row["ref_reads"]) / row["total_reads"]) * 100, 1)
            if row["total_reads"] > 0
            else 0,
            axis=1,
        )
        pileupfile["called_alt_perc"] = pileupfile.apply(
            lambda row: round((row["highest_alt_reads"] / row["total_reads"]) * 100, 1)
            if row["total_reads"] > 0
            else 0,
            axis=1,
        )

        pileupfile["called_base"] = pileupfile.apply(
            lambda row: row["refbase"]
            if row["called_ref_perc"] >= row["called_alt_perc"]
            else row["highest_alt_reads_base"],
            axis=1,
        )
        pileupfile["called_perc"] = pileupfile.apply(
            lambda row: row["called_ref_perc"]
            if row["called_ref_perc"] >= row["called_alt_perc"]
            else row["called_alt_perc"],
            axis=1,
        ).astype(float)
        pileupfile["called_reads"] = pileupfile.apply(
            lambda row: row["ref_reads"]
            if row["called_ref_perc"] >= row["called_alt_perc"]
            else row["highest_alt_reads"],
            axis=1,
        ).astype(int)

        intersect_pos = np.intersect1d(pileupfile["pos"], markerfile["pos"])
        markerfile = markerfile.loc[markerfile["pos"].isin(intersect_pos)]
        markerfile = markerfile.sort_values(by=["pos"])
        pileupfile = pileupfile.loc[pileupfile["pos"].isin(intersect_pos)]

        pileupfile = pileupfile.drop(["chr"], axis=1)
        df = pd.merge(markerfile, pileupfile, on="pos")

        df["state"] = df.apply(
            lambda row: "A"
            if row["called_base"] == row["anc"]
            else "D"
            if row["called_base"] == row["der"]
            else "NA",
            axis=1,
        )
        df["bool_state"] = df.apply(
            lambda row: True
            if (row["called_base"] == row["anc"] or row["called_base"] == row["der"])
            else False,
            axis=1,
        )

        markerfile_len = len(markerfile)

        # valid markers from positionsfile.txt
        general_info_list = ["Total of mapped reads: VCF", "Total of unmapped reads: VCF"]
        general_info_list += ["Valid markers: " + str(markerfile_len)]

        index_belowzero = df[df["called_reads"] == 0].index
        # Use copy() to avoid chained assignment warning
        df_belowzero = df[df.index.isin(index_belowzero)].copy()
        df_belowzero = df_belowzero.drop(["refbase", "altbase"], axis=1)
        df_belowzero["called_perc"] = "NA"
        df_belowzero["called_base"] = "NA"
        df_belowzero["state"] = "NA"
        df_belowzero["Description"] = "Position with zero reads"

        df = df[~df.index.isin(index_belowzero)]
        bool_list_state = df[df["bool_state"] == False].index

        # discordant genotypes
        # Use copy()
        df_discordantgenotype = df[df.index.isin(bool_list_state)].copy()
        df_discordantgenotype = df_discordantgenotype.drop(["bool_state"], axis=1)
        df_discordantgenotype["state"] = "NA"
        df_discordantgenotype["Description"] = "Discordant genotype"
        df = df[~df.index.isin(bool_list_state)]

        reads_thresh = int(args.reads_treshold)

        # read threshold
        # Use copy() because we filter df and then assign
        df_readsthreshold = df[df["called_reads"] < reads_thresh].copy()
        df_readsthreshold["Description"] = "Below read threshold"
        df = df[df["called_reads"] >= reads_thresh]

        base_majority = float(args.base_majority)

        # filter by base percentage
        df_basemajority = df[df["called_perc"] < base_majority].copy()
        df_basemajority["Description"] = "Below base majority"
        df = df[df["called_perc"] >= base_majority]

        df_fmf = pd.concat(
            [df_belowzero, df_readsthreshold, df_basemajority, df_discordantgenotype],
            axis=0,
            sort=True,
        )
        df_fmf["reads"] = df_fmf["called_reads"]
        df_fmf = df_fmf[
            [
                "chr",
                "pos",
                "marker_name",
                "haplogroup",
                "mutation",
                "anc",
                "der",
                "reads",
                "called_perc",
                "called_base",
                "state",
                "Description",
            ]
        ]

        df_out = df.copy()
        df_out["reads"] = df_out["called_reads"]
        df_out = df_out[
            [
                "chr",
                "pos",
                "marker_name",
                "haplogroup",
                "mutation",
                "anc",
                "der",
                "reads",
                "called_perc",
                "called_base",
                "state",
            ]
        ]

        general_info_list.append("Markers with zero reads: " + str(len(df_belowzero)))
        general_info_list.append(
            "Markers below the read threshold {"
            + str(reads_thresh)
            + "}: "
            + str(len(df_readsthreshold))
        )
        general_info_list.append(
            "Markers below the base majority threshold {"
            + str(base_majority)
            + "}: "
            + str(len(df_basemajority))
        )
        general_info_list.append(
            "Markers with discordant genotype: " + str(len(df_discordantgenotype))
        )
        general_info_list.append(
            "Markers without haplogroup information: " + str(len(df_fmf))
        )
        general_info_list.append("Markers with haplogroup information: " + str(len(df_out)))

        self.write_info_file(sample_vcf_folder, general_info_list)

        use_old = args.use_old
        outputfile = sample_vcf_folder / (sample_vcf_file.name.replace(".vcf.gz", ".out"))
        fmf_output = sample_vcf_folder / (sample_vcf_file.name.replace(".vcf.gz", ".fmf"))

        if use_old:
            df_out = df_out.sort_values(by=["haplogroup"], ascending=True)
            df_out = df_out[
                [
                    "chr",
                    "pos",
                    "marker_name",
                    "haplogroup",
                    "mutation",
                    "anc",
                    "der",
                    "reads",
                    "called_perc",
                    "called_base",
                    "state",
                ]
            ]
            df_fmf.to_csv(fmf_output, sep="\t", index=False)
            df_out.to_csv(outputfile, sep="\t", index=False)
            return

        df_out = df_out[
            [
                "chr",
                "pos",
                "marker_name",
                "haplogroup",
                "mutation",
                "anc",
                "der",
                "reads",
                "called_perc",
                "called_base",
                "state",
            ]
        ]
        df_fmf.to_csv(fmf_output, sep="\t", index=False)

        # sort based on the tree
        lst_df = df_out.values.tolist()
        mappable_df = {}
        for lst in lst_df:
            if lst[3] not in mappable_df:
                mappable_df[lst[3]] = []
            mappable_df[lst[3]].append(lst)

        tree = Tree(
            self.config.data_folder
            / "hg_prediction_tables"
            / yleaf_constants.TREE_FILE
        )
        with open(outputfile, "w") as f:
            f.write(
                "\t".join(
                    [
                        "chr",
                        "pos",
                        "marker_name",
                        "haplogroup",
                        "mutation",
                        "anc",
                        "der",
                        "reads",
                        "called_perc",
                        "called_base",
                        "state",
                        "depth\n",
                    ]
                )
            )
            for node_key in tree.node_mapping:
                if node_key not in mappable_df:
                    continue
                depth = tree.get(node_key).depth
                for lst in mappable_df[node_key]:
                    f.write("\t".join(map(str, lst)) + f"\t{depth}\n")

        LOG.info(f"Finished extracting genotypes for {sample_vcf_file.name}")

    def main_vcf_split(
        self,
        position_bed_file: Path,
        base_out_folder: Path,
        args: argparse.Namespace,
        vcf_file: Path,
    ):
        # first sort the vcf file
        sorted_vcf_file = base_out_folder / (
            vcf_file.name.replace(".vcf.gz", ".sorted.vcf.gz")
        )
        cmd = f"bcftools sort -O z -o {sorted_vcf_file} {vcf_file}"
        try:
            run_command(cmd, shell=True)
        except ExternalCommandError:
            LOG.error(f"Failed to sort the vcf file {vcf_file.name}. Skipping...")
            return None

        # next index the vcf file
        cmd = f"bcftools index -f {sorted_vcf_file}"
        run_command(cmd, shell=True)

        # get chromosome annotation
        cmd = f"bcftools query -f '%CHROM\n' {sorted_vcf_file} | uniq"
        # Must capture output
        result = run_command(cmd, shell=True, capture_output=True, text=True)
        chromosomes = result.stdout.strip().split("\n")
        chry = [x for x in chromosomes if "y" in x.lower()]

        if len(chry) == 0:
            LOG.error("Unable to find Y-chromosome in the vcf file. Exiting...")
            raise YleafError("Unable to find Y-chromosome in the vcf file.")
        elif len(chry) > 1:
            LOG.error("Multiple Y-chromosome annotations found in the vcf file. Exiting...")
            raise YleafError("Multiple Y-chromosome annotations found in the vcf file.")
        else:
            # make new position_bed_file with correct chrY annotation
            new_position_bed_file = base_out_folder / (
                vcf_file.name.replace(".vcf.gz", "temp_position_bed.bed")
            )
            with open(position_bed_file) as f:
                with open(new_position_bed_file, "w") as f2:
                    for line in f:
                        line = line.replace("chrY", chry[0])
                        f2.write(line)

        # filter the vcf file using the reference bed file
        filtered_vcf_file = (
            base_out_folder
            / "filtered_vcf_files"
            / (sorted_vcf_file.name.replace(".sorted.vcf.gz", ".filtered.vcf.gz"))
        )
        cmd = f"bcftools view -O z -R {new_position_bed_file} {sorted_vcf_file} > {filtered_vcf_file}"
        run_command(cmd, shell=True)

        # remover temp_position_bed.bed
        cmd = ["rm", str(new_position_bed_file)]
        run_command(cmd)

        # remove sorted.vcf.gz and sorted.vcf.gz.csi
        run_command(["rm", str(sorted_vcf_file)])
        run_command(["rm", f"{sorted_vcf_file}.csi"])

        # check number of samples in the vcf file
        cmd = f"bcftools query -l {filtered_vcf_file} | wc -l"
        result = run_command(cmd, shell=True, capture_output=True)
        num_samples = int(result.stdout.strip())

        if num_samples > 1:
            # split the vcf file into separate files for each sample
            split_vcf_folder = base_out_folder / (
                vcf_file.name.replace(".vcf.gz", "_split")
            )
            self.safe_create_dir(split_vcf_folder, args.force)
            cmd = f"bcftools +split {filtered_vcf_file} -Oz -o {split_vcf_folder}"
            run_command(cmd, shell=True)
            sample_vcf_files = self.get_files_with_extension(split_vcf_folder, ".vcf.gz")
        elif num_samples == 1:
            sample_vcf_files = [filtered_vcf_file]
        else:
            LOG.error("No samples found in the vcf file. Exiting...")
            raise YleafError("No samples found in the vcf file.")

        return sample_vcf_files

    def main_vcf(self, args: argparse.Namespace, base_out_folder: Path):
        position_bed_file = self.get_position_bed_file(
            args.reference_genome, args.use_old, args.ancient_DNA
        )
        path_markerfile = self.get_position_file(
            args.reference_genome, args.use_old, args.ancient_DNA
        )

        self.safe_create_dir(base_out_folder / "filtered_vcf_files", args.force)

        files = self.get_files_with_extension(args.vcffile, ".vcf.gz")

        if not args.reanalyze:
            with multiprocessing.Pool(processes=args.threads) as p:
                sample_vcf_files = p.map(
                    partial(self.main_vcf_split, position_bed_file, base_out_folder, args), files
                )

            sample_vcf_files = [x for x in sample_vcf_files if x is not None]
            sample_vcf_files = [item for sublist in sample_vcf_files for item in sublist]

            with multiprocessing.Pool(processes=args.threads) as p:
                p.map(
                    partial(self.run_vcf, path_markerfile, base_out_folder, args),
                    sample_vcf_files,
                )

        else:
            with multiprocessing.Pool(processes=args.threads) as p:
                p.map(partial(self.run_vcf, path_markerfile, base_out_folder, args), files)

    def main_fastq(self, args: argparse.Namespace, out_folder: Path):
        files = self.get_files_with_extension(args.fastq, ".fastq")
        files += self.get_files_with_extension(args.fastq, ".fastq.gz")
        reference = self.get_reference_path(args.reference_genome, True)
        bam_folder = out_folder / yleaf_constants.FASTQ_BAM_FILE_FOLDER
        try:
            os.mkdir(bam_folder)
        except OSError:
            pass
        LOG.info("Creating bam files from fastq files...")

        # For paired end fastq files (with _R1 and _R2 in the file name) and gzipped fastq files with .gz extension align the pairs together
        for fastq_file in files:
            if (
                "_R1" in str(fastq_file)
                and Path(str(fastq_file).replace("_R1", "_R2")).exists()
            ):
                fastq_file2 = Path(str(fastq_file).replace("_R1", "_R2"))
                LOG.info(f"Starting with running for {fastq_file} and {fastq_file2}")
                sam_file = bam_folder / "temp_fastq_sam.sam"
                fastq_cmd = f"minimap2 -ax sr -k14 -w7 -t {args.threads} {reference} {fastq_file} {fastq_file2} > {sam_file}"
                run_command(fastq_cmd, shell=True)
                bam_file = bam_folder / (fastq_file.name.rsplit("_R1", 1)[0] + ".bam")
                cmd = f"samtools view -@ {args.threads} -bS {sam_file} | samtools sort -@ {args.threads} -m 2G -o {bam_file}"
                run_command(cmd, shell=True)
                cmd = f"samtools index -@ {args.threads} {bam_file}"
                run_command(cmd, shell=True)
                os.remove(sam_file)
            elif (
                "_R1" not in str(fastq_file)
                and "_R2" not in str(fastq_file)
                and ".gz" in str(fastq_file)
            ):
                LOG.info(f"Starting with running for {fastq_file}")
                sam_file = bam_folder / "temp_fastq_sam.sam"

                fastq_cmd = f"minimap2 -ax sr -k14 -w7 -t {args.threads} {reference} {fastq_file} > {sam_file}"
                run_command(fastq_cmd, shell=True)
                bam_file = bam_folder / (fastq_file.name.rsplit(".", 1)[0] + ".bam")
                cmd = f"samtools view -@ {args.threads} -bS {sam_file} | samtools sort -@ {args.threads} -m 2G -o {bam_file}"
                run_command(cmd, shell=True)
                cmd = f"samtools index -@ {args.threads} {bam_file}"
                run_command(cmd, shell=True)
                os.remove(sam_file)
        args.bamfile = bam_folder
        self.main_bam_cram(args, out_folder, True)

    def main_bam_cram(self, args: argparse.Namespace, base_out_folder: Path, is_bam: bool):
        if args.bamfile is not None:
            files = self.get_files_with_extension(args.bamfile, ".bam")
        elif args.cramfile is not None:
            if args.cram_reference is None:
                raise ValueError("Please specify a reference genome for the CRAM file.")
            files = self.get_files_with_extension(args.cramfile, ".cram")
        else:
            print("Please specify either (a) bam or a cram file(s)")
            return

        with multiprocessing.Pool(processes=args.threads) as p:
            p.map(partial(self.run_bam_cram, args, base_out_folder, is_bam), files)

    def run_bam_cram(
        self,
        args: argparse.Namespace, base_out_folder: Path, is_bam: bool, input_file: Path
    ):
        LOG.info(f"Starting with running for {input_file}")
        output_dir = base_out_folder / input_file.name.rsplit(".", 1)[0]
        self.safe_create_dir(output_dir, args.force)
        general_info_list = self.samtools(output_dir, input_file, is_bam, args)
        self.write_info_file(output_dir, general_info_list)
        if args.private_mutations:
            self.find_private_mutations(output_dir, input_file, args, is_bam)
        LOG.debug(f"Finished running for {input_file.name}")
        print()

    def get_files_with_extension(self, path: str | Path, ext: str) -> list[Path]:
        """Get all files with a certain extension from a path. The path can be a file or a dir."""
        filtered_files = []
        path = Path(path)  # to be sure
        if path.is_dir():
            for file in path.iterdir():
                if str(file)[-len(ext) :] == ext:
                    filtered_files.append(file)
            return filtered_files
        else:
            return [path]

    def check_reference(self, requested_version: str):
        reference_file = self.get_reference_path(requested_version, True)
        if not reference_file.exists() or os.path.getsize(reference_file) < 100:
            LOG.info(
                f"No reference genome version was found. Downloading the {requested_version} reference genome. This "
                f"should be a one time thing."
            )
            download_reference.main(requested_version)
            LOG.info("Finished downloading the reference genome.")

    def get_reference_path(self, requested_version: str, is_full: bool) -> Path | None:
        if is_full:
            if requested_version == yleaf_constants.HG19:
                reference_file = self.config.hg19_full_genome
            else:
                reference_file = self.config.hg38_full_genome
        else:
            if requested_version == yleaf_constants.HG19:
                reference_file = self.config.hg19_y_chromosome
            else:
                reference_file = self.config.hg38_y_chromosome
        return reference_file

    def samtools(
        self,
        output_folder: Path,
        path_file: Path,
        is_bam_pathfile: bool,
        args: argparse.Namespace,
    ) -> list[str]:
        outputfile = output_folder / (output_folder.name + ".out")
        fmf_output = output_folder / (output_folder.name + ".fmf")
        pileupfile = output_folder / "temp_haplogroup_pileup.pu"
        reference = args.cram_reference if not is_bam_pathfile else None

        if is_bam_pathfile:
            if not any(
                [
                    Path(str(path_file) + ".bai").exists(),
                    Path(str(path_file).rstrip(".bam") + ".bai").exists(),
                ]
            ):
                cmd = f"samtools index -@{args.threads} {path_file}"
                run_command(cmd, shell=True)
        else:
            if not any(
                [
                    Path(str(path_file) + ".crai").exists(),
                    Path(str(path_file).rstrip(".cram") + ".crai").exists(),
                ]
            ):
                cmd = f"samtools index -@{args.threads} {path_file}"
                run_command(cmd, shell=True)
        header, mapped, unmapped = self.chromosome_table(
            path_file, output_folder, output_folder.name
        )

        position_file = self.get_position_file(
            args.reference_genome, args.use_old, args.ancient_DNA
        )

        bed = output_folder / "temp_position_bed.bed"
        self.write_bed_file(bed, position_file, header)

        self.execute_mpileup(bed, path_file, pileupfile, args.quality_thresh, reference)

        general_info_list = [
            "Total of mapped reads: " + str(mapped),
            "Total of unmapped reads: " + str(unmapped),
        ]

        self.extract_haplogroups(
            position_file,
            args.reads_treshold,
            args.base_majority,
            pileupfile,
            fmf_output,
            outputfile,
            is_bam_pathfile,
            args.use_old,
            general_info_list,
        )

        os.remove(pileupfile)
        os.remove(bed)

        LOG.debug("Finished extracting haplogroups")
        return general_info_list

    def chromosome_table(
        self, path_file: Path, path_folder: Path, file_name: str
    ) -> tuple[str, int, int]:
        output = path_folder / (file_name + ".chr")
        cmd = ["samtools", "idxstats", str(path_file)]

        LOG.debug(f"Started running the following command: {' '.join(cmd)}")
        # Streaming to pandas
        # samtools idxstats output is usually small, so communicate() is fine and avoids deadlocks.
        # But for consistency and robustness, let's stick to Popen/communicate which reads both buffers.
        try:
            with subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as process:
                stdout_data, stderr_data = process.communicate()
                if process.returncode != 0:
                     raise ExternalCommandError(" ".join(cmd), process.returncode, stderr_data.decode("utf-8"))
        except Exception as e:
            raise ExternalCommandError(" ".join(cmd), -1, str(e))

        df_chromosome = pd.read_csv(io.BytesIO(stdout_data), sep="\t", header=None)

        total_reads = sum(df_chromosome[2])

        # Check if Y exists
        y_rows = df_chromosome[df_chromosome[0].str.contains("Y")]
        if y_rows.empty:
            # Handle case where no Y
             unmapped = 0 # or error
        else:
             unmapped = y_rows[3].values[0]

        df_chromosome["perc"] = (df_chromosome[2] / total_reads) * 100
        df_chromosome = df_chromosome.round(decimals=2)
        df_chromosome["perc"] = df_chromosome["perc"].astype(str) + "%"
        df_chromosome = df_chromosome.drop(columns=[1, 3])
        df_chromosome.columns = ["chr", "reads", "perc"]
        df_chromosome.to_csv(output, index=None, sep="\t")

        if "Y" in df_chromosome["chr"].values:
            return "Y", total_reads, unmapped
        elif "chrY" in df_chromosome["chr"].values:
            return "chrY", total_reads, unmapped
        else:
            LOG.error("Unable to find Y-chromosomal data in the provided files. Exiting...")
            raise YleafError("No Y-chromosomal data")

    def get_position_file(
        self,
        reference_name: str,
        use_old: bool,
        ancient_DNA: bool,
    ) -> Path:
        if use_old:
            if ancient_DNA:
                position_file = (
                    self.config.data_folder
                    / reference_name
                    / yleaf_constants.OLD_POSITION_ANCIENT_FILE
                )
            else:
                position_file = (
                    self.config.data_folder
                    / reference_name
                    / yleaf_constants.OLD_POSITION_FILE
                )
        else:
            if ancient_DNA:
                position_file = (
                    self.config.data_folder
                    / reference_name
                    / yleaf_constants.NEW_POSITION_ANCIENT_FILE
                )
            else:
                position_file = (
                    self.config.data_folder
                    / reference_name
                    / yleaf_constants.NEW_POSITION_FILE
                )
        return position_file

    def get_position_bed_file(
        self,
        reference_name: str,
        use_old: bool,
        ancient_DNA: bool,
    ) -> Path:
        if use_old:
            if ancient_DNA:
                position_file = (
                    self.config.data_folder
                    / reference_name
                    / yleaf_constants.OLD_POSITION_ANCIENT_BED_FILE
                )
            else:
                position_file = (
                    self.config.data_folder
                    / reference_name
                    / yleaf_constants.OLD_POSITION_BED_FILE
                )
        else:
            if ancient_DNA:
                position_file = (
                    self.config.data_folder
                    / reference_name
                    / yleaf_constants.NEW_POSITION_ANCIENT_BED_FILE
                )
            else:
                position_file = (
                    self.config.data_folder
                    / reference_name
                    / yleaf_constants.NEW_POSITION_BED_FILE
                )
        return position_file

    def write_bed_file(self, bed: Path, markerfile: Path, header: str):
        mf = pd.read_csv(markerfile, sep="\t", header=None)
        mf = mf[[0, 3]]
        mf[0] = header
        mf.to_csv(str(bed), sep="\t", index=False, header=False)

    def execute_mpileup(
        self,
        bed: Path | None,
        bam_file: Path,
        pileupfile: Path,
        quality_thresh: float,
        reference: Path | None,
    ):
        cmd = "samtools mpileup"
        if bed is not None:
            cmd += f" -l {str(bed)}"

        if reference is not None:
            cmd += f" -f {str(reference)}"
        cmd += f" -AQ{quality_thresh}q1 {str(bam_file)} > {str(pileupfile)}"
        run_command(cmd, shell=True)

    def extract_haplogroups(
        self,
        path_markerfile: Path,
        reads_thresh: float,
        base_majority: int,
        path_pileupfile: Path,
        fmf_output: Path,
        outputfile: Path,
        is_bam_file: bool,
        use_old: bool,
        general_info_list: list[str],
    ):
        LOG.debug("Starting with extracting haplogroups...")
        markerfile = pd.read_csv(path_markerfile, header=None, sep="\t")
        markerfile.columns = [
            "chr",
            "marker_name",
            "haplogroup",
            "pos",
            "mutation",
            "anc",
            "der",
        ]
        markerfile = markerfile.drop_duplicates(subset="pos", keep="first", inplace=False)

        # packagemanagement is the best
        pileupfile = pd.read_csv(
            path_pileupfile,
            header=None,
            sep="\t",
            dtype={0: str, 1: int, 2: str, 3: int, 4: str, 5: str},
            on_bad_lines="skip",
        )

        pileupfile.columns = ["chr", "pos", "refbase", "reads", "align", "quality"]

        if not is_bam_file:
            ref_base = pileupfile["refbase"].values
            read_results = pileupfile["align"].values
            new_read_results = list(map(self.replace_with_bases, ref_base, read_results))
            pileupfile["align"] = new_read_results

        intersect_pos = np.intersect1d(pileupfile["pos"], markerfile["pos"])
        markerfile = markerfile.loc[markerfile["pos"].isin(intersect_pos)]
        markerfile = markerfile.sort_values(by=["pos"])
        pileupfile = pileupfile.loc[pileupfile["pos"].isin(intersect_pos)]

        pileupfile = pileupfile.drop(["chr"], axis=1)
        df = pd.merge(markerfile, pileupfile, on="pos")

        markerfile_len = len(markerfile)

        # valid markers from positionsfile.txt
        general_info_list.append("Valid markers: " + str(markerfile_len))

        index_belowzero = df[df["reads"] == 0].index
        # Use copy()
        df_belowzero = df[df.index.isin(index_belowzero)].copy()
        df_belowzero = df_belowzero.drop(["refbase", "align", "quality"], axis=1)
        df_belowzero["called_perc"] = "NA"
        df_belowzero["called_base"] = "NA"
        df_belowzero["state"] = "NA"
        df_belowzero["Description"] = "Position with zero reads"

        df = df[~df.index.isin(index_belowzero)]

        freq_dict = self.get_frequency_table(df.values)
        df_freq_table = pd.DataFrame.from_dict(freq_dict, orient="index")
        df_freq_table.columns = ["A", "T", "G", "C", "+", "-"]

        df_freq_table = df_freq_table.drop(["+", "-"], axis=1)
        df = df.drop(["refbase", "align", "quality"], axis=1)

        list_col_indices = np.argmax(df_freq_table.values, axis=1)
        called_base = df_freq_table.columns[list_col_indices]  # noqa
        total_count_bases = np.sum(df_freq_table.values, axis=1)
        max_count_bases = np.max(df_freq_table, axis=1)
        called_perc = round((max_count_bases / total_count_bases) * 100, 1)

        bool_anc = np.equal(np.array(called_base), df["anc"].values)
        bool_der = np.equal(np.array(called_base), df["der"].values)

        bool_list_anc = np.where(bool_anc, "A", "D")
        bool_list_anc = bool_list_anc.astype("object")
        bool_list_der = np.where(bool_der, "D", "A")
        bool_list_der = bool_list_der.astype("object")
        bool_list_state = np.equal(bool_list_anc, bool_list_der)

        df["called_perc"] = np.array(called_perc, dtype=int)
        df["called_base"] = called_base
        df["state"] = bool_list_anc
        df["bool_state"] = bool_list_state

        # discordant genotypes
        # Use copy()
        df_discordantgenotype = df[~bool_list_state].copy()
        df_discordantgenotype = df_discordantgenotype.drop(["bool_state"], axis=1)
        df_discordantgenotype["state"] = "NA"
        df_discordantgenotype["Description"] = "Discordant genotype"
        df = df[bool_list_state]

        # read threshold
        # Use copy()
        df_readsthreshold = df[df["reads"] < reads_thresh].copy()
        df_readsthreshold["Description"] = "Below read threshold"
        df = df[df["reads"] >= reads_thresh]

        # filter by base percentage
        # Use copy()
        df_basemajority = df[df["called_perc"] < base_majority].copy()
        df_basemajority["Description"] = "Below base majority"
        df = df[df["called_perc"] >= base_majority]

        df_fmf = pd.concat(
            [df_belowzero, df_readsthreshold, df_basemajority, df_discordantgenotype],
            axis=0,
            sort=True,
        )
        df_fmf = df_fmf[
            [
                "chr",
                "pos",
                "marker_name",
                "haplogroup",
                "mutation",
                "anc",
                "der",
                "reads",
                "called_perc",
                "called_base",
                "state",
                "Description",
            ]
        ]

        df_out = df.drop(["bool_state"], axis=1)

        general_info_list.append("Markers with zero reads: " + str(len(df_belowzero)))
        general_info_list.append(
            "Markers below the read threshold {"
            + str(reads_thresh)
            + "}: "
            + str(len(df_readsthreshold))
        )
        general_info_list.append(
            "Markers below the base majority threshold {"
            + str(base_majority)
            + "}: "
            + str(len(df_basemajority))
        )
        general_info_list.append(
            "Markers with discordant genotype: " + str(len(df_discordantgenotype))
        )
        general_info_list.append(
            "Markers without haplogroup information: " + str(len(df_fmf))
        )
        general_info_list.append("Markers with haplogroup information: " + str(len(df_out)))

        if use_old:
            df_out = df_out.sort_values(by=["haplogroup"], ascending=True)
            df_out = df_out[
                [
                    "chr",
                    "pos",
                    "marker_name",
                    "haplogroup",
                    "mutation",
                    "anc",
                    "der",
                    "reads",
                    "called_perc",
                    "called_base",
                    "state",
                ]
            ]
            df_fmf.to_csv(fmf_output, sep="\t", index=False)
            df_out.to_csv(outputfile, sep="\t", index=False)
            return

        df_out = df_out[
            [
                "chr",
                "pos",
                "marker_name",
                "haplogroup",
                "mutation",
                "anc",
                "der",
                "reads",
                "called_perc",
                "called_base",
                "state",
            ]
        ]
        df_fmf.to_csv(fmf_output, sep="\t", index=False)

        # sort based on the tree
        lst_df = df_out.values.tolist()
        mappable_df = {}
        for lst in lst_df:
            if lst[3] not in mappable_df:
                mappable_df[lst[3]] = []
            mappable_df[lst[3]].append(lst)

        tree = Tree(
            self.config.data_folder
            / "hg_prediction_tables"
            / yleaf_constants.TREE_FILE
        )
        with open(outputfile, "w") as f:
            f.write(
                "\t".join(
                    [
                        "chr",
                        "pos",
                        "marker_name",
                        "haplogroup",
                        "mutation",
                        "anc",
                        "der",
                        "reads",
                        "called_perc",
                        "called_base",
                        "state",
                        "depth\n",
                    ]
                )
            )
            for node_key in tree.node_mapping:
                if node_key not in mappable_df:
                    continue
                depth = tree.get(node_key).depth
                for lst in mappable_df[node_key]:
                    f.write("\t".join(map(str, lst)) + f"\t{depth}\n")


    def replace_with_bases(self, base: str, read_result: str) -> str:
        return read_result.replace(",", base[0]).replace(".", base[0])


    def get_frequency_table(self, mpileup: list[Any]) -> dict[str, list[int]]:
        frequency_table = {}
        for i in mpileup:
            fastadict = self.get_frequencies(i[9])
            frequency_table[i[3]] = list(fastadict.values())
        return frequency_table


    def get_frequencies(self, sequence: str) -> dict[str, int]:
        sequence = sequence.upper()

        # 1. Remove start markers "^."
        sequence = START_RE.sub("", sequence)

        # 2. Remove end markers "$"
        sequence = sequence.replace("$", "")

        # 3. Handle Indels
        indel_counts = {"+": 0, "-": 0}

        clean_parts = []
        last_pos = 0

        # We assume pileup indel sequences (AGCTN*) do not contain + or - followed by digits.
        for match in INDEL_RE.finditer(sequence):
            start, end = match.span()

            # Keep the part before the indel
            if start > last_pos:
                clean_parts.append(sequence[last_pos:start])

            # Count the indel
            indel_type = match.group(1)
            indel_counts[indel_type] += 1

            # Get length to skip
            skip_len = int(match.group(2))

            # Advance last_pos over the match AND the skipped sequence
            last_pos = end + skip_len

        # Append remaining part
        if last_pos < len(sequence):
            clean_parts.append(sequence[last_pos:])

        final_seq = "".join(clean_parts)

        # 4. Count bases
        # We only care about A, T, G, C, *, and we map * to -
        c = Counter(final_seq)

        result = {
            "A": c["A"],
            "T": c["T"],
            "G": c["G"],
            "C": c["C"],
            "-": c["*"] + indel_counts["-"],
            "+": indel_counts["+"],
        }

        return result


    def write_info_file(self, folder: Path, general_info_list: list[str]):
        try:
            with open(folder / (folder.name + ".info"), "a") as f:
                for marker in general_info_list:
                    f.write(marker)
                    f.write("\n")
        except OSError:
            LOG.warning("Failed to write .info file")


    def find_private_mutations(
        self, output_folder: Path, path_file: Path, args: argparse.Namespace, is_bam: bool
    ):
        # identify mutations not part of haplogroups that are annotated in dbsnp or differ from the reference genome
        LOG.debug("Starting with extracting private mutations...")
        snp_reference_file = self.get_reference_path(args.reference_genome, False)
        snp_database_file = (
            self.config.data_folder
            / args.reference_genome
            / yleaf_constants.SNP_DATA_FILE
        )

        # run mpileup
        pileup_file = output_folder / "temp_private_mutation_pileup.pu"
        self.execute_mpileup(
            None,
            path_file,
            pileup_file,
            args.quality_thresh,
            snp_reference_file if not is_bam else None,
        )

        LOG.debug("Loading reference files")

        position_file = self.get_position_file(
            args.reference_genome, args.use_old, args.ancient_DNA
        )
        filter_positions = self.load_filter_data(position_file)
        known_snps = self.load_snp_database_file(snp_database_file, args.minor_allele_frequency)
        ychrom_reference = self.load_reference_file(snp_reference_file)

        LOG.debug("Finding private mutations...")
        private_mutations = []
        confirmed_private_mutations = []
        with open(pileup_file) as f:
            for index, line in enumerate(f):
                try:
                    chrom, position, ref_base, count, aligned, quality = (
                        line.strip().split()
                    )
                except ValueError:
                    LOG.warning(f"failed to read line {index} of pileupfile")
                    continue
                if chrom != "chrY":
                    continue
                # not enough reads
                if int(count) < args.reads_treshold:
                    continue
                if not is_bam:
                    aligned = self.replace_with_bases(ref_base, aligned)
                if position in filter_positions:
                    continue

                # not annotated in dbsnp
                if position not in known_snps and snp_reference_file is not None:
                    freq_dict = self.get_frequencies(aligned)
                    actual_allele, allele_count = max(freq_dict.items(), key=lambda x: x[1])
                    called_percentage = round(
                        allele_count / sum(freq_dict.values()) * 100, 2
                    )
                    # not a high enough base majority measured, cannot be sure of real allele
                    if called_percentage < args.base_majority:
                        continue
                    ref_base = ychrom_reference[int(position) - 1]

                    # do not match against insertions or repeat regions (lower case)
                    if (
                        ref_base not in ACCEPTED_REF_BASES
                        or actual_allele not in ACCEPTED_REF_BASES
                    ):
                        continue
                    if ref_base == actual_allele:
                        continue
                    private_mutations.append(
                        f"{chrom}\t{position}\t-\t{ref_base}->{actual_allele}\t{ref_base}\t"
                        f"{actual_allele}\t{allele_count}\t{called_percentage}\tNA\n"
                    )
                elif position in known_snps and snp_database_file is not None:
                    freq_dict = self.get_frequencies(aligned)
                    actual_allele, allele_count = max(freq_dict.items(), key=lambda x: x[1])
                    called_percentage = round(
                        allele_count / sum(freq_dict.values()) * 100, 2
                    )
                    # not a high enough base majority measured, cannot be sure of real allele
                    if called_percentage < args.base_majority:
                        continue
                    possible_minor_alleles = {
                        dct["minor_allele"] for dct in known_snps[position]
                    }
                    # measured allele is not the dbsnp allele
                    if actual_allele not in possible_minor_alleles:
                        continue

                    matched_pos_dct = [
                        dct
                        for dct in known_snps[position]
                        if dct["minor_allele"] == actual_allele
                    ][0]
                    rs, major_allele, minor_allele, frequency = matched_pos_dct.values()
                    confirmed_private_mutations.append(
                        f"{chrom}\t{position}\t{rs}\t{major_allele}->{actual_allele}\t"
                        f"{major_allele}\t{actual_allele}\t{allele_count}\t"
                        f"{called_percentage}\t{frequency}\n"
                    )
        os.remove(pileup_file)
        with open(output_folder / f"{output_folder.name}.pmu", "w") as f:
            f.write(
                "chrom\tposition\trn_no\tmutation\treference\tdetected\treads\tcalled_percentage\t"
                "minor allele frequency\n"
            )
            f.write("".join(confirmed_private_mutations))
            f.write("".join(private_mutations))
        LOG.debug("Finished extracting private mutations")


    def load_filter_data(self, path: Path) -> set[str]:
        if self.cached_position_data is None:
            self.cached_position_data = set()
            with open(path) as f:
                for line in f:
                    self.cached_position_data.add(line.strip().split("\t")[3])
        return self.cached_position_data


    def load_snp_database_file(
        self, path: Path, minor_allele_frequency: float
    ) -> dict[str, list[dict[str, str]]] | None:
        if path is not None:
            if self.cached_snp_database is None:
                self.cached_snp_database = defaultdict(list)
                with open(path) as f:
                    f.readline()
                    for line in f:
                        rs, position, major_allele, minor_allele, frequency = (
                            line.strip().split(",")
                        )
                        if float(frequency) > minor_allele_frequency:
                            continue
                        self.cached_snp_database[position].append(
                            {
                                "rs": rs,
                                "major_allele": major_allele,
                                "minor_allele": minor_allele,
                                "frequency": frequency,
                            }
                        )
                self.cached_snp_database = dict(self.cached_snp_database)
            return self.cached_snp_database
        else:
            return None


    def load_reference_file(self, path: Path) -> str | None:
        if path is not None:
            if self.cached_reference_file is None:
                parts = []
                with open(path) as f:
                    for line in f:
                        if line.startswith(">"):
                            continue
                        parts.append(line.strip())
                self.cached_reference_file = "".join(parts)
            return self.cached_reference_file
        else:
            return None


    def run_haplogroup_prediction(
        self,
        path_file: Path,
        output: Path,
        use_old: bool,
        prediction_quality: float,
        threads: int,
    ):
        if use_old:
            script = yleaf_constants.SRC_FOLDER / "old_predict_haplogroup.py"
            cmd = f"python {script} -i {path_file} -o {output}"
            run_command(cmd, shell=True)
        else:
            from yleaf import predict_haplogroup

            namespace = argparse.Namespace(
                input=path_file,
                outfile=output,
                minimum_score=prediction_quality,
                threads=threads,
            )
            predict_haplogroup.main(namespace)


    def draw_haplogroups(self, haplogroup_file: Path, collapsed_draw_mode: bool):
        # make sure that it is only imported if requested by user
        from yleaf import haplogroup_tree_image

        namespace = argparse.Namespace(
            input=haplogroup_file,
            collapse_mode=collapsed_draw_mode,
            outfile=haplogroup_file.parent / HAPLOGROUP_IMAGE_FILE_NAME,
        )
        haplogroup_tree_image.main(namespace)

    def safe_create_dir(self, folder: Path, force: bool):
        """Create the given folder. If the folder is already present delete if the user agrees."""
        if folder.is_dir():
            while True and not force:
                LOG.warning(
                    "Folder "
                    + str(folder)
                    + " already exists, would you like to remove it?"
                )
                choice = input("y/n: ")
                if str(choice).upper() == "Y":
                    break
                elif str(choice).upper() == "N":
                    sys.exit(0) # Should we raise exception instead? Main logic typically handles user exit.
                    # But replacing sys.exit with exception is better for testing.
                    # However, this is interactive input.
                else:
                    print("Please type y/Y or n/N")
            shutil.rmtree(folder)
            os.mkdir(folder)
        else:
            try:
                os.mkdir(folder)
            except OSError:
                print("Failed to create directory. Exiting...")
                raise YleafError("Failed to create directory")


def setup_logger(out_folder: Path):
    """Setup logging"""
    LOG.setLevel(logging.DEBUG)

    # Remove existing handlers to avoid duplicates if called multiple times
    for handler in LOG.handlers[:]:
        LOG.removeHandler(handler)

    handler = logging.StreamHandler(sys.stdout)
    start_time = time.time()
    formatter = MyFormatter(
        "%(levelname)s %(currentTime)s (%(passedTime)s s) - %(message)s",
        starttime=start_time,
    )
    handler.setFormatter(formatter)
    handler.setLevel(logging.INFO)
    LOG.addHandler(handler)

    file_handler = logging.FileHandler(filename=out_folder / "run.log")
    file_handler.setFormatter(formatter)
    file_handler.setLevel(logging.DEBUG)
    LOG.addHandler(file_handler)

    LOG.debug("Logger created")


def get_arguments() -> argparse.Namespace:
    parser = ArgumentParser()

    parser.add_argument(
        "-fastq",
        "--fastq",
        required=False,
        help="Use raw FastQ files",
        metavar="PATH",
        type=check_file,
    )
    parser.add_argument(
        "-bam",
        "--bamfile",
        required=False,
        help="input BAM file",
        metavar="PATH",
        type=check_file,
    )
    parser.add_argument(
        "-cram",
        "--cramfile",
        required=False,
        help="input CRAM file",
        metavar="PATH",
        type=check_file,
    )
    parser.add_argument(
        "-cr",
        "--cram_reference",
        required=False,
        help="Reference genome for the CRAM file. Required when using CRAM files.",
        metavar="PATH",
        type=check_file,
    )
    parser.add_argument(
        "-vcf",
        "--vcffile",
        required=False,
        help="input VCF file (.vcf.gz)",
        metavar="PATH",
        type=check_file,
    )
    parser.add_argument(
        "-ra",
        "--reanalyze",
        required=False,
        help="reanalyze (skip filtering and splitting) the vcf file",
        action="store_true",
    )
    parser.add_argument(
        "-force", "--force", action="store_true", help="Delete files without asking"
    )
    parser.add_argument(
        "-rg",
        "--reference_genome",
        help="The reference genome build to be used. If no reference is available "
        "they will be downloaded. If you added references in your config.txt file these"
        " will be used instead as reference or the location will be used to download the "
        "reference if those files are missing or empty.",
        choices=[yleaf_constants.HG19, yleaf_constants.HG38],
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Folder name containing outputs",
        metavar="STRING",
    )
    parser.add_argument(
        "-r",
        "--reads_treshold",
        help="The minimum number of reads for each base. (default=10)",
        type=int,
        required=False,
        default=10,
        metavar="INT",
    )
    parser.add_argument(
        "-q",
        "--quality_thresh",
        help="Minimum quality for each read, integer between 10 and 40. [10-40] (default=20)",
        type=int,
        required=False,
        metavar="INT",
        default=20,
    )
    parser.add_argument(
        "-b",
        "--base_majority",
        help="The minimum percentage of a base result for acceptance, integer between 50 and 99."
        " [50-99] (default=90)",
        type=int,
        required=False,
        metavar="INT",
        default=90,
    )
    parser.add_argument(
        "-t",
        "--threads",
        dest="threads",
        help="The number of processes to use when running Yleaf.",
        type=int,
        default=1,
        metavar="INT",
    )
    parser.add_argument(
        "-pq",
        "--prediction_quality",
        type=float,
        required=False,
        default=0.95,
        metavar="FLOAT",
        help="The minimum quality of the prediction (QC-scores) for it to be accepted. [0-1] (default=0.95)",
    )

    # arguments for prediction
    parser.add_argument(
        "-old",
        "--use_old",
        dest="use_old",
        help="Add this value if you want to use the old prediction method of Yleaf (version 2.3). This"
        " version only uses the ISOGG tree and slightly different prediction criteria.",
        action="store_true",
    )

    # arguments for drawing haplo group trees
    parser.add_argument(
        "-dh",
        "--draw_haplogroups",
        help="Draw the predicted haplogroups in the haplogroup tree.",
        action="store_true",
    )
    parser.add_argument(
        "-hc",
        "--collapsed_draw_mode",
        help="Add this flag to compress the haplogroup tree image and remove all uninformative "
        "haplogroups from it.",
        action="store_true",
    )

    # arguments for ancient DNA samples
    parser.add_argument(
        "-aDNA",
        "--ancient_DNA",
        help="Add this flag if the sample is ancient DNA. This will ignore all G > A and C > T mutations.",
        action="store_true",
    )

    # arguments for private mutations
    parser.add_argument(
        "-p",
        "--private_mutations",
        help="Add this flag to search for private mutations. These are variations that are not"
        " considered in the phylogenetic tree and thus not used for haplogroup prediction, "
        "however can be informative and differentiate individuals within the same haplogroup "
        "prediction.",
        action="store_true",
    )

    parser.add_argument(
        "-maf",
        "--minor_allele_frequency",
        help="Maximum rate of minor allele for it to be considered"
        " as a private mutation. (default=0.01)",
        default=0.01,
        type=float,
        metavar="FLOAT",
    )

    args = parser.parse_args()
    return args


def check_file(path: str) -> Path:
    """Check for the presence of a file and return a Path object"""
    object_path = Path(path)
    if not object_path.exists():
        raise argparse.ArgumentTypeError("Path to provided file/dir does not exist")
    return object_path


def logo():
    print(r"""

                   |
                  /|\
                 /\|/\
                \\\|///
                 \\|//
                  |||
                  |||
                  |||

        """)


def main():
    print(
        "Erasmus MC Department of Genetic Identification\nYleaf: software tool for human Y-chromosomal "
        f"phylogenetic analysis and haplogroup inference v{__version__}"
    )
    logo()

    try:
        args = get_arguments()
        out_folder = Path(args.output)

        config = Configuration()
        pipeline = YleafPipeline(config)

        # Create output folder (recreates if exists and user agrees/force)
        pipeline.safe_create_dir(out_folder, args.force)

        setup_logger(out_folder)
        LOG.info(f"Running Yleaf with command: {' '.join(sys.argv)}")

        # make sure the reference genome is present before doing something else, if not present it is downloaded
        pipeline.check_reference(args.reference_genome)

        if args.fastq:
            pipeline.main_fastq(args, out_folder)
        elif args.bamfile:
            pipeline.main_bam_cram(args, out_folder, True)
        elif args.cramfile:
            pipeline.main_bam_cram(args, out_folder, False)
        elif args.vcffile:
            pipeline.main_vcf(args, out_folder)
        else:
            LOG.error("Please specify either a bam, a cram, a fastq, or a vcf file")
            raise ValueError("Please specify either a bam, a cram, a fastq, or a vcf file")

        hg_out = out_folder / PREDICTION_OUT_FILE_NAME
        pipeline.run_haplogroup_prediction(
            out_folder, hg_out, args.use_old, args.prediction_quality, args.threads
        )
        if args.draw_haplogroups:
            pipeline.draw_haplogroups(hg_out, args.collapsed_draw_mode)

        LOG.info("Done!")

    except YleafError as e:
        LOG.error(f"Yleaf error: {e}")
        sys.exit(1)
    except Exception as e:
        LOG.error(f"Unexpected error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
