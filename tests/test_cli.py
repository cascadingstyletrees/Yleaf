import itertools
import sys

import pytest

from yleaf import Yleaf


@pytest.fixture
def populated_dir(tmp_path):
    folder = tmp_path / "output"
    folder.mkdir()
    (folder / "keep.txt").write_text("keep")
    return folder


@pytest.mark.parametrize(
    ("flag", "attribute"),
    [
        ("-fastq", "fastq"),
        ("--fastq", "fastq"),
        ("-bam", "bamfile"),
        ("--bamfile", "bamfile"),
        ("-cram", "cramfile"),
        ("--cramfile", "cramfile"),
        ("-vcf", "vcffile"),
        ("--vcffile", "vcffile"),
    ],
)
def test_each_input_flag_is_accepted(flag, attribute, tmp_path, monkeypatch):
    sample = tmp_path / "sample"
    sample.touch()
    monkeypatch.setattr(
        sys, "argv", ["Yleaf", flag, str(sample), "-o", "output", "-rg", "hg38"]
    )

    args = Yleaf.get_arguments()

    assert getattr(args, attribute) == sample
    assert args.reference_genome == "hg38"
    assert args.threads == 1
    assert args.force is False


def test_missing_input_is_rejected_before_output_handling(populated_dir, monkeypatch):
    monkeypatch.setattr(
        sys, "argv", ["Yleaf", "-o", str(populated_dir), "-rg", "hg38"]
    )
    monkeypatch.setattr("builtins.input", lambda _: pytest.fail("unexpected prompt"))
    monkeypatch.setattr(
        Yleaf, "check_reference", lambda _: pytest.fail("unexpected reference setup")
    )

    with pytest.raises(SystemExit) as error:
        Yleaf.main()

    assert error.value.code == 2
    assert (populated_dir / "keep.txt").read_text() == "keep"


@pytest.mark.parametrize(
    ("first", "second"),
    list(itertools.combinations(["-bam", "-cram", "-fastq", "-vcf"], 2)),
)
def test_multiple_input_types_are_rejected(first, second, tmp_path, monkeypatch):
    sample = tmp_path / "sample"
    sample.touch()
    monkeypatch.setattr(
        sys,
        "argv",
        ["Yleaf", first, str(sample), second, str(sample), "-o", "out", "-rg", "hg38"],
    )

    with pytest.raises(SystemExit) as error:
        Yleaf.get_arguments()

    assert error.value.code == 2


@pytest.mark.parametrize("choice", ["", "n", "N", "no", " NO "])
def test_declining_preserves_output(choice, populated_dir, monkeypatch):
    responses = iter([choice])
    monkeypatch.setattr("builtins.input", lambda _: next(responses))

    with pytest.raises(SystemExit) as error:
        Yleaf.safe_create_dir(populated_dir, False)

    assert error.value.code == 0
    assert (populated_dir / "keep.txt").read_text() == "keep"


@pytest.mark.parametrize("choice", ["y", "Y", "yes", " YES "])
def test_confirmation_recreates_output(choice, populated_dir, monkeypatch):
    responses = iter([choice])
    monkeypatch.setattr("builtins.input", lambda _: next(responses))

    Yleaf.safe_create_dir(populated_dir, False)

    assert populated_dir.is_dir()
    assert list(populated_dir.iterdir()) == []


def test_force_recreates_output_without_prompt(populated_dir, monkeypatch):
    monkeypatch.setattr("builtins.input", lambda _: pytest.fail("unexpected prompt"))

    Yleaf.safe_create_dir(populated_dir, True)

    assert populated_dir.is_dir()
    assert list(populated_dir.iterdir()) == []


def test_new_output_is_created_without_prompt(tmp_path, monkeypatch):
    folder = tmp_path / "new-output"
    monkeypatch.setattr("builtins.input", lambda _: pytest.fail("unexpected prompt"))

    Yleaf.safe_create_dir(folder, False)

    assert folder.is_dir()
