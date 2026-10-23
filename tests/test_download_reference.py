import gzip
import io

import pytest

from yleaf import download_reference


@pytest.fixture
def progress(monkeypatch):
    now = [0.0]
    monkeypatch.setattr(
        download_reference, "monotonic", lambda: now[0], raising=False
    )
    return now


def test_progress_is_throttled_and_final_update_is_not_delayed(progress, capsys):
    reporter = download_reference.DownloadProgress()
    reporter(0, 8192, 8192000)
    for block in range(1, 100):
        progress[0] = block / 1000
        reporter(block, 8192, 8192000)
    progress[0] = 1.1
    reporter(100, 8192, 8192000)
    reporter(1000, 8192, 8192000)

    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err.count("Downloading:") == 3
    assert "100.0%" in captured.err
    assert "\r" not in captured.err


def test_progress_clamps_overshoot_and_does_not_repeat_completion(progress, capsys):
    reporter = download_reference.DownloadProgress()
    reporter(2, 8192, 10000)
    reporter(1, 10000, 10000)

    captured = capsys.readouterr()
    assert captured.err.count("100.0%") == 1
    assert "163.8%" not in captured.err


def test_unknown_size_still_reports_downloaded_bytes(progress, capsys):
    reporter = download_reference.DownloadProgress()
    reporter(0, 1024 * 1024, -1)
    progress[0] = 1.1
    reporter(2, 1024 * 1024, -1)

    captured = capsys.readouterr()
    assert "2.0 MiB" in captured.err
    assert "%" not in captured.err


def test_terminal_progress_rewrites_the_line(progress, monkeypatch):
    terminal = io.StringIO()
    monkeypatch.setattr(terminal, "isatty", lambda: True)
    monkeypatch.setattr(download_reference.sys, "stderr", terminal)

    download_reference.DownloadProgress()(0, 8192, 10000)

    assert terminal.getvalue().startswith("\rDownloading:")
    assert not terminal.getvalue().endswith("\n")


def test_download_hook_and_reference_extraction(tmp_path, monkeypatch, capsys):
    reference = tmp_path / "reference.fa"
    chromosome = tmp_path / "chrY.fa"
    reference.touch()
    chromosome.touch()
    monkeypatch.setattr(
        download_reference.yleaf_constants, "HG19_FULL_GENOME", reference
    )
    monkeypatch.setattr(
        download_reference.yleaf_constants, "HG19_Y_CHROMOSOME", chromosome
    )
    data = ">chr1\n" + "A" * 120 + "\n>chrY\n" + "C" * 120 + "\n>chr2\nGG\n"

    def retrieve(url, destination, reporthook):
        assert url.endswith("/hg19.fa.gz")
        reporthook(0, 8192, -1)
        with gzip.open(destination, "wb") as output:
            output.write(data.encode())
        reporthook(1, 8192, -1)
        return str(destination), {}

    monkeypatch.setattr(download_reference.urllib.request, "urlretrieve", retrieve)

    download_reference.install_genome_files("hg19")

    assert reference.read_text() == data
    assert chromosome.read_text() == ">chrY\n" + "C" * 120 + "\n"
    assert not reference.with_suffix(".fa.gz").exists()
    captured = capsys.readouterr()
    assert captured.out == ""
    assert "100.0%" in captured.err
