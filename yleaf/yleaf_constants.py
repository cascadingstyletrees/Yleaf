from pathlib import Path

SRC_FOLDER: Path = Path(__file__).absolute().parent
DATA_FOLDER: Path = SRC_FOLDER / "data"
# CONFIG_PATH is handled by Configuration class now
CHR_NAMING_CONVENTION_FILE: Path = SRC_FOLDER / "chr_naming_convention.txt"
HG_PREDICTION_FOLDER: Path = DATA_FOLDER / "hg_prediction_tables"

HG19: str = "hg19"
HG38: str = "hg38"

FULL_REF_FILE: str = "full_reference.fa"
Y_REF_FILE: str = "chrY.fa"
SNP_DATA_FILE: str = "snp_data.csv"
NEW_POSITION_FILE: str = "new_positions.txt"
OLD_POSITION_FILE: str = "old_positions.txt"
NEW_POSITION_BED_FILE: str = "new_positions.bed"
OLD_POSITION_BED_FILE: str = "old_positions.bed"
NEW_POSITION_ANCIENT_FILE: str = "new_positions_ancient.txt"
OLD_POSITION_ANCIENT_FILE: str = "old_positions_ancient.txt"
NEW_POSITION_ANCIENT_BED_FILE: str = "new_positions_ancient.bed"
OLD_POSITION_ANCIENT_BED_FILE: str = "old_positions_ancient.bed"

TREE_FILE: str = "tree.json"

FASTQ_BAM_FILE_FOLDER: str = "bam_files"
