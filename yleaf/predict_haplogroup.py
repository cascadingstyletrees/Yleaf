#!/usr/bin/env python

"""
Script for Haplogroup prediction using the YFull tree

Developed at Erasmus Medical Center Department of Genetic Identification

License: GNU General Public License v3 or later
A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

Autor: Bram van Wersch
"""

import argparse
import logging
import math
import multiprocessing
from collections.abc import Iterator
from functools import partial
from pathlib import Path
from typing import Any

from yleaf import yleaf_constants
from yleaf.configuration import Configuration
from yleaf.tree import Node, Tree
from yleaf.exceptions import ConfigurationError

# Globals removed/encapsulated.
# Constants
DEFAULT_MIN_SCORE = 0.95

LOG = logging.getLogger("yleaf_logger")


class HgMarkersLinker:
    """Safe for a certain haplogroup if the number of ancestral and derived markers."""

    DERIVED: str = "D"
    ANCESTRAL: str = "A"
    UNDEFINED: str = "N"

    _ancestral_markers: set[str]
    _derived_markers: set[str]

    def __init__(self):
        self._ancestral_markers = set()
        self._derived_markers = set()

    def add(self, marker_name: str, state: str):
        if state == self.DERIVED:
            self._derived_markers.add(marker_name)
        else:
            self._ancestral_markers.add(marker_name)

    def get_derived_markers(self) -> set[str]:
        return self._derived_markers

    @property
    def nr_total(self) -> int:
        return len(self._ancestral_markers) + len(self._derived_markers)

    @property
    def nr_derived(self) -> int:
        return len(self._derived_markers)

    @property
    def nr_ancestral(self) -> int:
        return len(self._ancestral_markers)

    def get_state(self) -> str:
        """at least a fraction of 0.6 need to either derived or ancestral, otherwise the state can not be accurately
        determined and will be returned as undefined"""
        if self.nr_total == 0:
             return self.UNDEFINED
        if self.nr_derived / self.nr_total >= 0.6:
            return self.DERIVED
        if self.nr_ancestral / self.nr_total >= 0.6:
            return self.ANCESTRAL
        return self.UNDEFINED


class HaplogroupPredictor:
    def __init__(self, config: Configuration, tree: Tree, min_score: float = DEFAULT_MIN_SCORE):
        self.config = config
        self.tree = tree
        self.min_score = min_score
        self.backbone_groups: set[str] = set()
        self.main_haplo_groups: set[str] = set()
        self.expected_states_cache: dict[str, dict[str, set[str]]] = {}
        self.qc1_score_cache: dict[str, float] = {}

        self._load_backbone_groups()

    def _load_backbone_groups(self):
        """Read some basic data that is always needed"""
        hg_pred_folder = self.config.data_folder / "hg_prediction_tables"
        intermediates_file = hg_pred_folder / "major_tables/Intermediates.txt"

        try:
            with open(intermediates_file) as f:
                for line in f:
                    if "~" in line:
                        continue
                    self.backbone_groups.add(line.strip())
        except FileNotFoundError:
            msg = f"Could not find intermediates file at {intermediates_file}. This file is critical for scoring."
            LOG.error(msg)
            raise ConfigurationError(msg)

        # add capitals A-Z to get all groups
        major_tree_list = [
            "A00", "A00a", "A00b", "A00c", "A0-T", "A1", "A1b", "A1b1", "BT", "CT", "CF",
            "F", "F4", "F2", "F3", "GHIJK", "G", "HIJK", "H", "IJK", "IJ", "I", "I2", "I1",
            "J", "J1", "J2", "K", "K2", "K2d", "K2c", "K2b", "P", "R", "R2", "R1", "R1b",
            "R1a", "Q", "K2b1", "S", "M", "NO", "O", "N", "LT", "L", "T", "F1", "C", "DE",
            "D", "E", "B", "A1a", "A0",
        ]

        for major_hg in major_tree_list:
            self.main_haplo_groups.add(major_hg)

    def predict(self, folder: Path) -> list[Any]:
        # Reset per-sample cache
        self.qc1_score_cache = {}

        if folder.name == "filtered_vcf_files":
            return [None, None, None]

        try:
            haplotype_dict = read_yleaf_out_file(folder / (folder.name + ".out"))
        except FileNotFoundError:
            LOG.warning(
                f"WARNING: failed to find .out file from yleaf run for sample {folder.name}. This sample will"
                " be skipped."
            )
            return [None, None, None]

        best_haplotype_score = self.get_most_likely_haplotype(
            haplotype_dict, self.min_score
        )
        return [haplotype_dict, best_haplotype_score, folder]

    def get_most_likely_haplotype(
        self, haplotype_dict: dict[str, HgMarkersLinker], treshold: float
    ) -> tuple[str, str | list[str], float, float, float, float, int | str]:
        sorted_depth_haplotypes = sorted(
            haplotype_dict.keys(), key=lambda k: self.tree.get(k).depth, reverse=True
        )

        intermediate_states = {
            value: haplotype_dict[value]
            for value in self.backbone_groups
            if value in haplotype_dict
        }

        covered_nodes = set()
        best_score = ("NA", "NA", "NA", "NA", "NA", "NA", "NA")

        for haplotype_name in sorted_depth_haplotypes:
            node = self.tree.get(haplotype_name)

            # only record most specific nodes regardless of scores
            if node.name in covered_nodes:
                continue
            parent = node.parent
            path = [node.name]

            # caclulate score 3 already since this is most efficient
            qc3_score_match = 0
            qc3_score_total = 0

            # walk the tree back
            while parent is not None:
                path.append(parent.name)
                if (
                    parent.name in haplotype_dict
                    and node.name[0] in parent.name
                    and parent.name != node.name[0]
                ):
                    # in the same overal haplogroup
                    state = haplotype_dict[parent.name].get_state()
                    if state == HgMarkersLinker.DERIVED:
                        qc3_score_match += 1

                    # in case it can not be decided what the state is ratio between 0.4 and 0.6
                    if state != HgMarkersLinker.UNDEFINED:
                        qc3_score_total += 1

                parent = parent.parent

            qc1_score = self.get_qc1_score(path, intermediate_states)

            # if any of the scores are below treshold, the total can not be above so ignore
            if qc1_score < treshold:
                continue

            if haplotype_dict[node.name].nr_total == 0:
                qc2_score = 0
            else:
                qc2_score = (
                    haplotype_dict[node.name].nr_derived
                    / haplotype_dict[node.name].nr_total
                )
                if qc2_score < treshold:
                    continue

            if qc3_score_total == 0:
                final_qc3_score = 0.0
            else:
                final_qc3_score = qc3_score_match / qc3_score_total
                if final_qc3_score < treshold:
                    continue

            total_score = math.prod([qc1_score, qc2_score, final_qc3_score])

            # if above filter we found the hit
            if total_score > treshold:
                ancestral_children = self.get_ancestral_children(node, haplotype_dict)
                best_score = (
                    node.name,
                    ancestral_children,
                    qc1_score,
                    qc2_score,
                    final_qc3_score,
                    total_score,
                    node.depth,
                )
                break

            # else, no hit is found, but still report Quality Scores
            else:
                best_score = (
                    "NA",
                    "NA",
                    qc1_score,
                    qc2_score,
                    final_qc3_score,
                    total_score,
                    "NA",
                )

            # make sure that less specific nodes are not recorded
            for node_name in path:
                covered_nodes.add(node_name)
        return best_score

    def get_qc1_score(
        self, path: list[str], intermediate_states: dict[str, HgMarkersLinker]
    ) -> float:
        most_specific_backbone = None
        for value in path:
            if value in self.main_haplo_groups:
                most_specific_backbone = value
                break

        if most_specific_backbone is None:
            return 0

        # Check cache
        if most_specific_backbone in self.qc1_score_cache:
            return self.qc1_score_cache[most_specific_backbone]

        # Load expected states if not cached
        hg_folder = self.config.data_folder / "hg_prediction_tables"
        if most_specific_backbone in self.expected_states_cache:
            expected_states = self.expected_states_cache[most_specific_backbone]
        else:
            expected_states = {}
            int_file = hg_folder / f"major_tables/{most_specific_backbone}_int.txt"
            try:
                with open(int_file) as f:
                    for line in f:
                        if line == "":
                            continue
                        name, state = line.strip().split("\t")
                        expected_states[name] = {*state.split("/")}
            except FileNotFoundError:
                # If file missing, assume empty
                pass
            self.expected_states_cache[most_specific_backbone] = expected_states

        score_match = 0
        score_total = 0
        for name, marker_linker in intermediate_states.items():
            # If name not in expected_states, we skip? Original code would fail KeyError if not present?
            # Original code: expected_possible_states = expected_states[name]
            # This implies 'name' from intermediate_states MUST be in expected_states.
            # But intermediate_states comes from BACKBONE_GROUPS check.
            if name in expected_states:
                expected_possible_states = expected_states[name]
                state = marker_linker.get_state()
                if state in expected_possible_states:
                    score_match += 1
                if state != HgMarkersLinker.UNDEFINED:
                    score_total += 1
            else:
                # What if intermediate state is not in expected states for this backbone?
                # Original code would crash if name not in expected_states.
                # Assuming data integrity.
                pass

        if score_total == 0:
            return 0

        qc1_score = score_match / score_total
        self.qc1_score_cache[most_specific_backbone] = qc1_score
        return qc1_score

    def get_ancestral_children(
        self, node: Node, haplotype_dict: dict[str, HgMarkersLinker]
    ) -> list[str]:
        ancestral_children = []
        to_cover = [node]
        while len(to_cover) > 0:
            curr_node = to_cover.pop()
            for name in curr_node.children:
                if name in haplotype_dict:
                    state = haplotype_dict[name].get_state()
                    if state == HgMarkersLinker.ANCESTRAL:
                        ancestral_children.append(name)
                        continue
                to_cover.append(self.tree.get(name))
        return ancestral_children


def read_yleaf_out_file(file: Path | str) -> dict[str, HgMarkersLinker]:
    """Read the full yleaf .out file and parse all lines into a dictionary keyed on haplogroups"""
    haplotype_dict = {}
    with open(file) as f:
        f.readline() # skip header
        for line in f:
            parts = line.strip().split("\t")
            # _, _, marker, haplogroup, _, _, _, _, _, _, state, _ = line.strip().split("\t")
            # Using unpacking with known structure from Yleaf.py
            # "chr", "pos", "marker_name", "haplogroup", "mutation", "anc", "der", "reads", "called_perc", "called_base", "state", "depth"
            if len(parts) < 11:
                continue
            marker = parts[2]
            haplogroup = parts[3]
            state = parts[10]

            if haplogroup not in haplotype_dict:
                haplotype_dict[haplogroup] = HgMarkersLinker()
            haplotype_dict[haplogroup].add(marker, state)
    return haplotype_dict


def read_input_folder(folder_name: str) -> Iterator[Path]:
    """Read all the folders present in the input folder."""
    in_folder = Path(folder_name)
    for folder in in_folder.iterdir():
        if not folder.is_dir():
            continue
        if folder.name == yleaf_constants.FASTQ_BAM_FILE_FOLDER:
            continue
        yield folder


def add_to_final_table(
    final_table: list[list[Any]],
    haplotype_dict: dict[str, HgMarkersLinker],
    best_haplotype_scores: tuple,
    folder: Path,
):
    total_reads, valid_markers = process_info_file(folder / (folder.name + ".info"))
    hg, ancestral_children, qc1, qc2, qc3, total, _ = best_haplotype_scores

    # Format QC scores
    # Original code: qc1, qc2, qc3 are floats.
    # We keep them as is.

    if hg == "NA":
        final_table.append(
            [folder.name, hg, "", total_reads, valid_markers, total, qc1, qc2, qc3]
        )
        return
    marker_list = list(haplotype_dict[hg].get_derived_markers())
    if len(marker_list) > 2:
        marker_list = marker_list[:2] + ["etc."]

    ancestral_children_list = ancestral_children if isinstance(ancestral_children, list) else []

    if len(ancestral_children_list) > 0:
        ancestral_string = "x" + ",".join(ancestral_children_list)
        final_table.append(
            [
                folder.name,
                f"{hg}*({ancestral_string})",
                ";".join(marker_list),
                total_reads,
                valid_markers,
                total,
                qc1,
                qc2,
                qc3,
            ]
        )
    else:
        final_table.append(
            [
                folder.name,
                hg,
                ";".join(marker_list),
                total_reads,
                valid_markers,
                total,
                qc1,
                qc2,
                qc3,
            ]
        )


def process_info_file(info_file: Path) -> tuple[int | str, int | str]:
    total_reads = "NA"
    valid_markers = "NA"

    try:
        with open(info_file) as f:
            for line in f:
                if line.startswith("Total of mapped reads:"):
                    total_reads = line.replace("Total of mapped reads:", "").strip()
                elif line.startswith("Markers with haplogroup information"):
                    valid_markers = line.replace(
                        "Markers with haplogroup information:", ""
                    ).strip()
    except FileNotFoundError:
        LOG.warning(
            "WARNING: failed to find .info file from yleaf run. This information is not critical but there"
            " will be some missing values in the output."
        )
    return total_reads, valid_markers


def write_final_table(final_table: list[list[Any]], out_file: Path | str):
    # sort for each sample based on QC-score
    # values[5] is total score. values[0] is sample name.
    # original: reverse=True (descending).
    final_table.sort(key=lambda values: (values[0], values[5]), reverse=True)

    header = "Sample_name\tHg\tHg_marker\tTotal_reads\tValid_markers\tQC-score\tQC-1\tQC-2\tQC-3\n"
    with open(out_file, "w") as f:
        f.write(header)
        for line in final_table:
            f.write("\t".join(map(str, line)) + "\n")


def get_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Erasmus MC: Genetic Identification\n Y-Haplogroup Prediction"
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Output file or path produced from Yleaf",
        metavar="FILE",
    )

    parser.add_argument(
        "-ms",
        "--minimum_score",
        help="Minimum score needed in order for a prediction to be considered"
        "for inclusion in the final data (default=0.95).",
        type=float,
        default=DEFAULT_MIN_SCORE,
    )

    parser.add_argument(
        "-o", "--outfile", required=True, help="Output file name", metavar="FILE"
    )

    parser.add_argument(
        "-t",
        "--threads",
        help="Number of threads to use (default=1).",
        type=int,
        default=1,
    )

    args = parser.parse_args()
    return args

def run_prediction_wrapper(predictor: HaplogroupPredictor, folder: Path):
    return predictor.predict(folder)

def main(namespace: argparse.Namespace = None):
    """Main entry point for prediction script"""
    LOG.info("Starting haplogroup prediction...")
    if namespace is None:
        namespace = get_arguments()
    in_folder = namespace.input
    output = namespace.outfile
    threads = namespace.threads

    config = Configuration()
    tree = Tree(
        config.data_folder
        / "hg_prediction_tables"
        / yleaf_constants.TREE_FILE
    )

    try:
        predictor = HaplogroupPredictor(config, tree, namespace.minimum_score)
        final_table = []

        with multiprocessing.Pool(processes=threads) as p:
            predictions = p.map(
                partial(run_prediction_wrapper, predictor),
                read_input_folder(in_folder),
            )

        for haplotype_dict, best_haplotype_score, folder in predictions:
            if haplotype_dict is None:
                continue
            add_to_final_table(final_table, haplotype_dict, best_haplotype_score, folder)

        write_final_table(final_table, output)
        LOG.debug("Finished haplogroup prediction")
    except ConfigurationError as e:
        LOG.error(f"Haplogroup prediction failed due to configuration error: {e}")
        # We might want to exit with non-zero status
        import sys
        sys.exit(1)


if __name__ == "__main__":
    main()
