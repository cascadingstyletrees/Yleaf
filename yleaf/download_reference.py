#!/usr/bin/python
"""
Code for downloading the reference genome and extracting the specific y-chromomse data

Developed at Erasmus Medical Center Department of Genetic Identification

License: GNU General Public License v3 or later
A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

Autor: Diego Montiel Gonzalez
Extensively modified by: Bram van Wersch
"""

import gzip
import logging
import os
import shutil
import urllib.request
from pathlib import Path

from yleaf import yleaf_constants
from yleaf.configuration import Configuration

LOG: logging = logging.getLogger("yleaf_logger")


def main(choice: str):
    if choice == yleaf_constants.HG19:
        reference_choice = [yleaf_constants.HG19]
    elif choice == yleaf_constants.HG38:
        reference_choice = [yleaf_constants.HG38]
    else:
        reference_choice = [yleaf_constants.HG19, yleaf_constants.HG38]

    config = Configuration()
    # running downloader
    for dir_name in reference_choice:
        install_genome_files(dir_name, config)


def install_genome_files(reference_choice: str, config: Configuration):
    LOG.info(f"Starting with preparing {reference_choice}...")

    if reference_choice == yleaf_constants.HG19:
        ref_file = config.hg19_full_genome
    else:
        ref_file = config.hg38_full_genome

    # Ensure parent directory exists
    if not ref_file.parent.exists():
        ref_file.parent.mkdir(parents=True, exist_ok=True)

    ref_gz_file = Path(str(ref_file) + ".gz")
    try:
        # Check if file exists and is not empty (size > 100 bytes is arbitrary check from original)
        if (not ref_file.exists() or os.path.getsize(ref_file) < 100) and not ref_gz_file.exists():
            LOG.debug(f"Downloading the {reference_choice} genome...")
            urllib.request.urlretrieve(
                f"http://hgdownload.cse.ucsc.edu/goldenPath/{reference_choice}"
                f"/bigZips/{reference_choice}.fa.gz",
                ref_gz_file,
            )

        # Unpack if gz exists and ref file is missing/empty
        if ref_gz_file.exists() and (not ref_file.exists() or os.path.getsize(ref_file) < 100):
            LOG.debug("Unpacking the downloaded archive...")
            with gzip.open(ref_gz_file, "rb") as f_in:
                with open(ref_file, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(ref_gz_file)

        if reference_choice == yleaf_constants.HG19:
            ychrom_file = config.hg19_y_chromosome
        else:
            ychrom_file = config.hg38_y_chromosome

        if not ychrom_file.exists() or os.path.getsize(ychrom_file) < 100:
            LOG.debug("Writing Ychromosomal data")
            get_ychrom_data(ref_file, ychrom_file)
    # try and cleanup when user aborts, this attempts to not leave half downloaded files
    except KeyboardInterrupt:
        try:
            if ref_gz_file.exists():
                os.remove(ref_gz_file)
            if ref_file.exists():
                os.remove(ref_file)
        # skip on IOerrors and such
        finally:
            raise


def get_ychrom_data(full_data_path: Path, yhcrom_file: Path):
    with open(yhcrom_file, "w") as fo, open(full_data_path) as fi:
        record = False
        for line in fi:
            if line == ">chrY\n":
                record = True
                fo.write(line)
            elif record:
                if line.startswith(">"):
                    break
                fo.write(line)
