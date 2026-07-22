import gzip
import shutil
import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).parents[1]
SCRIPTS = ROOT / "scripts"


def run_script(name: str, *args: str, input_text: str = "") -> subprocess.CompletedProcess:
    return subprocess.run(
        [str(SCRIPTS / name), *args],
        input=input_text,
        text=True,
        capture_output=True,
        check=True,
    )


def test_scripts_have_valid_bash_syntax():
    scripts = sorted(str(path) for path in SCRIPTS.glob("*.sh"))
    subprocess.run(["bash", "-n", *scripts], check=True)


def test_fasta_len_counts_complete_sequences():
    result = run_script(
        "fasta_len.sh",
        "-",
        input_text=">seq1 description\nACGU\nAC\n>seq2\nuu\n",
    )

    assert result.stdout == "Length\tCount\n2\t1\n6\t1\n"


def test_fasta_u2t_preserves_header_and_converts_both_cases():
    result = run_script(
        "fasta_u2t.sh",
        "-",
        input_text=">seq1 description\nACGU\nuu\n",
    )

    assert result.stdout == ">seq1 description\nACGT\ntt\n"


def test_fastq_len_and_incomplete_record_error():
    result = run_script(
        "fastq_len.sh",
        "-",
        input_text="@r1\nACGT\n+\nIIII\n@r2\nAC\n+\nII\n",
    )
    assert result.stdout == "Length\tCount\n2\t1\n4\t1\n"

    invalid = subprocess.run(
        [str(SCRIPTS / "fastq_len.sh"), "-"],
        input="@r1\nACGT\n+\n",
        text=True,
        capture_output=True,
    )
    assert invalid.returncode != 0
    assert "incomplete four-line FASTQ record" in invalid.stderr


def test_get_chrom_length_reuses_existing_fai(tmp_path):
    reference = tmp_path / "reference.fa"
    reference.write_text(">chr1\nACGT\n>chr2\nAAA\n")
    Path(f"{reference}.fai").write_text("chr1\t4\t6\t4\t5\nchr2\t3\t17\t3\t4\n")
    output = tmp_path / "chrom.sizes"

    run_script("get_chrom_length.sh", str(reference), str(output))

    assert output.read_text() == "chr1\t4\nchr2\t3\n"


@pytest.mark.skipif(shutil.which("pigz") is None, reason="pigz is not installed")
def test_remove_duplicates_uses_paired_sequences(tmp_path):
    read1 = tmp_path / "R1.fq.gz"
    read2 = tmp_path / "R2.fq.gz"
    out1 = tmp_path / "out1.fq.gz"
    out2 = tmp_path / "out2.fq.gz"
    with gzip.open(read1, "wt") as handle:
        handle.write("@a\nACGT\n+\nIIII\n@b\nACGT\n+\nIIII\n@c\nAAAA\n+\nIIII\n")
    with gzip.open(read2, "wt") as handle:
        handle.write("@a\nTGCA\n+\nIIII\n@b\nTGCA\n+\nIIII\n@c\nCCCC\n+\nIIII\n")

    result = run_script(
        "remove_duplicates.sh",
        "2",
        str(read1),
        str(read2),
        str(out1),
        str(out2),
    )

    with gzip.open(out1, "rt") as handle:
        assert handle.read() == "@c\nAAAA\n+\nIIII\n@a\nACGT\n+\nIIII\n"
    with gzip.open(out2, "rt") as handle:
        assert handle.read() == "@c\nCCCC\n+\nIIII\n@a\nTGCA\n+\nIIII\n"
    assert "Input read pairs: 3" in result.stderr
    assert "Output read pairs: 2" in result.stderr
    assert "Duplication rate: 33.33%" in result.stderr
    assert not list(tmp_path.glob("out1.fq.gz.temp.*"))
