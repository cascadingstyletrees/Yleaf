import json
import logging
from pathlib import Path
from typing import Any

from yleaf.exceptions import ConfigurationError

LOG = logging.getLogger("yleaf_logger")

class Configuration:
    """
    Manages configuration for Yleaf, including file paths for reference genomes.
    """

    def __init__(self, config_path: Path | None = None):
        self._src_folder = Path(__file__).absolute().parent
        self._data_folder = self._src_folder / "data"
        self._config_path = config_path or (self._src_folder / "config.json")

        # Defaults
        hg19_folder = self._data_folder / "hg19"
        hg38_folder = self._data_folder / "hg38"

        self.hg19_full_genome = hg19_folder / "full_reference.fa"
        self.hg19_y_chromosome = hg19_folder / "chrY.fa"
        self.hg38_full_genome = hg38_folder / "full_reference.fa"
        self.hg38_y_chromosome = hg38_folder / "chrY.fa"

        self.load()

    def load(self):
        """Loads configuration from the file if it exists."""
        if not self._config_path.exists():
            # Check for legacy config.txt
            legacy_path = self._src_folder / "config.txt"
            if legacy_path.exists():
                LOG.info(f"Loading legacy configuration from {legacy_path}")
                self._load_legacy_config(legacy_path)
            return

        try:
            with open(self._config_path, "r") as f:
                data = json.load(f)
                self._apply_config(data)
        except json.JSONDecodeError as e:
            LOG.warning(f"Failed to parse config file {self._config_path}: {e}")
        except Exception as e:
            LOG.error(f"Error loading configuration: {e}")

    def _load_legacy_config(self, path: Path):
        """Parses the legacy key=value config.txt."""
        config_data = {}
        try:
            with open(path, "r") as f:
                for line in f:
                    if "=" not in line:
                        continue
                    name, value = line.strip().split("=", 1)
                    name = name.strip()
                    value = value.strip()
                    if not value:
                        continue

                    if name == "full hg19 genome fasta location":
                        config_data["hg19_full_genome"] = value
                    elif name == "full hg38 genome fasta location":
                        config_data["hg38_full_genome"] = value
                    elif name == "hg19 chromosome Y fasta location":
                        config_data["hg19_y_chromosome"] = value
                    elif name == "hg38 chromosome Y fasta location":
                        config_data["hg38_y_chromosome"] = value
            self._apply_config(config_data)
        except Exception as e:
            LOG.error(f"Error loading legacy configuration: {e}")

    def _apply_config(self, data: dict[str, Any]):
        """Applies configuration data to the instance."""
        if "hg19_full_genome" in data:
            self.hg19_full_genome = self._validate_path("hg19_full_genome", data["hg19_full_genome"])
        if "hg19_y_chromosome" in data:
            self.hg19_y_chromosome = self._validate_path("hg19_y_chromosome", data["hg19_y_chromosome"])
        if "hg38_full_genome" in data:
            self.hg38_full_genome = self._validate_path("hg38_full_genome", data["hg38_full_genome"])
        if "hg38_y_chromosome" in data:
            self.hg38_y_chromosome = self._validate_path("hg38_y_chromosome", data["hg38_y_chromosome"])

    def _validate_path(self, name: str, value: str) -> Path:
        """Validates and returns a Path object."""
        path = Path(value)
        # Note: The original code created the file if it didn't exist.
        # We will not do that automatically here to avoid side effects during load,
        # but we can check existence.
        # However, for now, we just return the Path. Creation logic should be explicit where needed.
        if path.suffix not in [".fa", ".fasta", ".fna"]:
             # This was a ValueError in original code
             LOG.warning(f"Configured path for {name} ({path}) does not match expected extensions (.fa, .fasta, .fna)")
        return path

    def save(self):
        """Saves the current configuration to the JSON file."""
        data = {
            "hg19_full_genome": str(self.hg19_full_genome),
            "hg19_y_chromosome": str(self.hg19_y_chromosome),
            "hg38_full_genome": str(self.hg38_full_genome),
            "hg38_y_chromosome": str(self.hg38_y_chromosome),
        }
        try:
            with open(self._config_path, "w") as f:
                json.dump(data, f, indent=4)
        except Exception as e:
            raise ConfigurationError(f"Failed to save configuration: {e}")

    @property
    def data_folder(self) -> Path:
        return self._data_folder
