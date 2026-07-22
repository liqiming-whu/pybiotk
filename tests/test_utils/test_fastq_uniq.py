import gzip

import pytest

from pybiotk.io import FastqPair
from pybiotk.utils import fastq_uniq
from pybiotk.utils.fastq_uniq import _make_key, main


def write_fastq(path, records):
    with open(path, "w") as fq:
        for name, sequence in records:
            fq.write(f"@{name}\n{sequence}\n+\n{'F' * len(sequence)}\n")


def read_fastq_names(path):
    with gzip.open(path, "rt") as fq:
        return [line[1:].strip() for index, line in enumerate(fq) if index % 4 == 0]


def test_single_end_filters_before_deduplicating_by_id(tmp_path):
    input_fq = tmp_path / "input.fq"
    output_fq = tmp_path / "output.fq.gz"
    write_fastq(input_fq, [("same", "AAA"), ("same", "ACGT"), ("other", "ACGT")])

    main(str(input_fq), str(output_fq), min_len=4, by="id")

    assert read_fastq_names(output_fq) == ["same", "other"]


def test_empty_input_is_supported(tmp_path):
    input_fq = tmp_path / "empty.fq"
    output_fq = tmp_path / "output.fq.gz"
    input_fq.touch()

    main(str(input_fq), str(output_fq))

    assert read_fastq_names(output_fq) == []


def test_summary_counts_short_and_duplicate_reads_separately(tmp_path, monkeypatch):
    input_fq = tmp_path / "input.fq"
    output_fq = tmp_path / "output.fq.gz"
    messages = []
    monkeypatch.setattr(fastq_uniq.logger, "info", messages.append)
    write_fastq(input_fq, [("short", "AAA"), ("first", "ACGT"), ("copy", "ACGT")])

    main(str(input_fq), str(output_fq), min_len=4)

    summary = messages[-2]
    assert read_fastq_names(output_fq) == ["first"]
    assert "reads too short (<4nt): 1" in summary
    assert "duplicate reads: 1" in summary
    assert "eligible reads: 2" in summary
    assert "duplication rate: 50.00%" in summary


def test_digest_preserves_field_boundaries():
    assert _make_key(("AC", "G"), "digest") != _make_key(("A", "CG"), "digest")
    assert len(_make_key(("A" * 150, "C" * 150), "digest")) == 16


def test_exact_key_mode_remains_available(tmp_path):
    input_fq = tmp_path / "input.fq"
    output_fq = tmp_path / "output.fq.gz"
    write_fastq(input_fq, [("first", "ACGT"), ("copy", "ACGT")])

    main(str(input_fq), str(output_fq), min_len=4, key_mode="exact")

    assert read_fastq_names(output_fq) == ["first"]


def test_gzip_compression_level_is_configurable(tmp_path):
    input_fq = tmp_path / "input.fq"
    output_fq = tmp_path / "output.fq.gz"
    write_fastq(input_fq, [("first", "ACGT")])

    main(str(input_fq), str(output_fq), min_len=4, compresslevel=1)

    assert output_fq.read_bytes()[8] == 4


def test_doc_contains_read_pair_collision_probabilities():
    assert "Unique read pairs  Key count" in fastq_uniq.__doc__
    assert "500,000,000        500,000,000   3.6734e-22" in fastq_uniq.__doc__
    assert "1,000,000,000      1,000,000,000 1.4694e-21" in fastq_uniq.__doc__
    assert "6.7534e-03 (0.6753%)" in fastq_uniq.__doc__
    assert "2.6741e-02 (2.6741%)" in fastq_uniq.__doc__
    assert "64-bit hash       about 74" in fastq_uniq.__doc__
    assert "saving            about 397 (81.4%)" in fastq_uniq.__doc__


@pytest.mark.parametrize(
    ("inputs", "outputs", "message"),
    [
        (["a", "b", "c"], ["x", "y", "z"], "one or two input"),
        (["a", "b"], ["x"], "number of output"),
    ],
)
def test_main_rejects_invalid_file_counts(inputs, outputs, message):
    with pytest.raises(ValueError, match=message):
        main(inputs, outputs)


def test_main_rejects_invalid_options():
    with pytest.raises(ValueError, match="min_len"):
        main("input.fq", "output.fq.gz", min_len=-1)
    with pytest.raises(ValueError, match="uniqueness key"):
        main("input.fq", "output.fq.gz", by="unknown")
    with pytest.raises(ValueError, match="key mode"):
        main("input.fq", "output.fq.gz", key_mode="unknown")
    with pytest.raises(ValueError, match="compresslevel"):
        main("input.fq", "output.fq.gz", compresslevel=10)


def test_pair_end_rejects_different_file_lengths(tmp_path):
    r1 = tmp_path / "r1.fq"
    r2 = tmp_path / "r2.fq"
    write_fastq(r1, [("a/1", "ACGT"), ("b/1", "TGCA")])
    write_fastq(r2, [("a/2", "ACGT")])

    with pytest.raises(RuntimeError, match="different numbers"):
        list(FastqPair(str(r1), str(r2)))


def test_pair_end_rejects_unmatched_names(tmp_path):
    r1 = tmp_path / "r1.fq"
    r2 = tmp_path / "r2.fq"
    write_fastq(r1, [("a/1", "ACGT")])
    write_fastq(r2, [("b/2", "ACGT")])

    with pytest.raises(RuntimeError, match="unmatched read names"):
        list(FastqPair(str(r1), str(r2)))
