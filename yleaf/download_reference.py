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
import sys
import urllib.request
from pathlib import Path
from time import monotonic

from yleaf import yleaf_constants

LOG: logging = logging.getLogger("yleaf_logger")


def main(choice: str):
    if choice == yleaf_constants.HG19:
        reference_choice = [yleaf_constants.HG19]
    elif choice == yleaf_constants.HG38:
        reference_choice = [yleaf_constants.HG38]
    else:
        reference_choice = [yleaf_constants.HG19, yleaf_constants.HG38]

    # running downloader
    for dir_name in reference_choice:
        install_genome_files(dir_name)


class DownloadProgress:
    """Report download progress at most once per second, plus completion."""

    def __init__(self):
        self.last_update = float("-inf")
        self.last_report = None

    def __call__(self, block_num, block_size, total_size):
        downloaded = max(0, block_num * block_size)
        complete = total_size > 0 and downloaded >= total_size
        if total_size > 0:
            downloaded = min(downloaded, total_size)
        report = (downloaded, total_size)
        now = monotonic()
        if report == self.last_report:
            return
        if not complete and now - self.last_update < 1:
            return
        self.last_update = now
        self.last_report = report

        detail = f"{downloaded / (1024 * 1024):.1f} MiB"
        if total_size > 0:
            detail = f"{downloaded * 100 / total_size:.1f}% ({detail})"
        terminal = sys.stderr.isatty()
        prefix = "\r" if terminal else ""
        end = "" if terminal else "\n"
        sys.stderr.write(f"{prefix}Downloading: {detail}{end}")
        sys.stderr.flush()


def install_genome_files(reference_choice: str):
    LOG.info(f"Starting with preparing {reference_choice}...")

    if reference_choice == yleaf_constants.HG19:
        ref_file = yleaf_constants.HG19_FULL_GENOME
    else:
        ref_file = yleaf_constants.HG38_FULL_GENOME

    ref_gz_file = Path(str(ref_file) + ".gz")
    try:
        if os.path.getsize(ref_file) < 100 and not ref_gz_file.exists():
            LOG.info(f"Downloading the {reference_choice} genome...")
            progress = DownloadProgress()
            try:
                urllib.request.urlretrieve(
                    f"http://hgdownload.cse.ucsc.edu/goldenPath/{reference_choice}"
                    f"/bigZips/{reference_choice}.fa.gz",
                    ref_gz_file,
                    reporthook=progress,
                )
                downloaded_size = ref_gz_file.stat().st_size
                progress(1, downloaded_size, downloaded_size)
            finally:
                if sys.stderr.isatty():
                    print(file=sys.stderr)
        if os.path.getsize(ref_file) < 100:
            LOG.debug("Unpacking the downloaded archive...")
            with gzip.open(ref_gz_file, "rb") as f_in:
                with open(ref_file, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(ref_gz_file)

        if reference_choice == yleaf_constants.HG19:
            ychrom_file = yleaf_constants.HG19_Y_CHROMOSOME
        else:
            ychrom_file = yleaf_constants.HG38_Y_CHROMOSOME

        if os.path.getsize(ychrom_file) < 100:
            LOG.debug("Writing Ychromosomal data")
            get_ychrom_data(ref_file, ychrom_file)
    # try and cleanup when user aborts, this attempts to not leave half downloaded files
    except KeyboardInterrupt:
        try:
            os.remove(ref_gz_file)
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
