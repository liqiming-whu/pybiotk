import gzip
import inspect

import pytest

from pybiotk.io import OpenFqGzip
from pybiotk.utils.fastx_rename import _format_index, fastx_rename, fastx_rename_pair


def write_fastq(path, records):
    with open(path, "w") as fq:
        for name, sequence in records:
            fq.write(f"@{name}\n{sequence}\n+\n{'F' * len(sequence)}\n")


def read_fastq_names(path):
    with gzip.open(path, "rt") as fq:
        return [line[1:].strip() for index, line in enumerate(fq) if index % 4 == 0]


def test_index_format_bases():
    assert _format_index(35, 36) == "z"
    assert _format_index(36, 36) == "10"
    assert _format_index(255, 16) == "ff"
    assert _format_index(9, 10) == "9"
    assert _format_index(10, 10) == "10"


def test_single_end_default_uses_short_decimal_names(tmp_path):
    input_fq = tmp_path / "input.fq"
    output_fq = tmp_path / "output.fq.gz"
    write_fastq(input_fq, [("same", "ACGT"), ("same", "TGCA")])

    fastx_rename(str(input_fq), str(output_fq))

    assert read_fastq_names(output_fq) == ["read_1", "read_2"]


def test_single_end_can_use_original_name(tmp_path):
    input_fq = tmp_path / "input.fq"
    output_fq = tmp_path / "output.fq.gz"
    write_fastq(input_fq, [("first", "ACGT"), ("second", "TGCA")])

    fastx_rename(str(input_fq), str(output_fq), use_original_name=True)

    assert read_fastq_names(output_fq) == ["first_1", "second_2"]


def test_single_end_custom_prefix(tmp_path):
    input_fq = tmp_path / "input.fq"
    output_fq = tmp_path / "output.fq.gz"
    write_fastq(input_fq, [("same", "ACGT"), ("same", "TGCA")])

    fastx_rename(str(input_fq), str(output_fq), prefix="r", index_base=16)

    assert read_fastq_names(output_fq) == ["r_1", "r_2"]


def test_single_end_preserve_mode(tmp_path):
    input_fq = tmp_path / "input.fq"
    output_fq = tmp_path / "output.fq.gz"
    write_fastq(input_fq, [("same", "ACGT"), ("same", "TGCA")])

    fastx_rename(str(input_fq), str(output_fq), mode="preserve")

    assert read_fastq_names(output_fq) == ["same", "same_2"]


def test_single_end_preserve_avoids_generated_suffix_collisions(tmp_path):
    input_fq = tmp_path / "input.fq"
    output_fq = tmp_path / "output.fq.gz"
    write_fastq(input_fq, [
        ("readA", "AAAA"),
        ("readA", "CCCC"),
        ("readA_2", "GGGG"),
        ("readA", "TTTT"),
    ])

    fastx_rename(str(input_fq), str(output_fq), mode="preserve")

    assert read_fastq_names(output_fq) == ["readA", "readA_2", "readA_2_2", "readA_3"]


def test_pair_index_default_uses_short_decimal_names(tmp_path):
    r1, r2 = tmp_path / "r1.fq", tmp_path / "r2.fq"
    out1, out2 = tmp_path / "out1.fq.gz", tmp_path / "out2.fq.gz"
    write_fastq(r1, [("orig/1", "ACGT"), ("other/1", "AAAA")])
    write_fastq(r2, [("orig/2", "TGCA"), ("other/2", "TTTT")])

    fastx_rename_pair([str(r1)], [str(r2)], str(out1), str(out2))

    assert read_fastq_names(out1) == ["read_1", "read_2"]
    assert read_fastq_names(out2) == ["read_1", "read_2"]


def test_pair_index_can_use_original_name_without_mate_suffix(tmp_path):
    r1, r2 = tmp_path / "r1.fq", tmp_path / "r2.fq"
    out1, out2 = tmp_path / "out1.fq.gz", tmp_path / "out2.fq.gz"
    write_fastq(r1, [("orig/1", "ACGT"), ("other/1", "AAAA")])
    write_fastq(r2, [("orig/2", "TGCA"), ("other/2", "TTTT")])

    fastx_rename_pair([str(r1)], [str(r2)], str(out1), str(out2),
                      use_original_name=True)

    assert read_fastq_names(out1) == ["orig_1", "other_2"]
    assert read_fastq_names(out2) == ["orig_1", "other_2"]


def test_pair_index_mode_keeps_global_index_across_files(tmp_path):
    r1a, r2a = tmp_path / "r1a.fq", tmp_path / "r2a.fq"
    r1b, r2b = tmp_path / "r1b.fq", tmp_path / "r2b.fq"
    out1, out2 = tmp_path / "out1.fq.gz", tmp_path / "out2.fq.gz"
    write_fastq(r1a, [("a/1", "ACGT")])
    write_fastq(r2a, [("a/2", "TGCA")])
    write_fastq(r1b, [("a/1", "AAAA"), ("b/1", "CCCC")])
    write_fastq(r2b, [("a/2", "TTTT"), ("b/2", "GGGG")])

    fastx_rename_pair([str(r1a), str(r1b)], [str(r2a), str(r2b)], str(out1), str(out2),
                      mode="index", prefix="pair", index_base=36)

    assert read_fastq_names(out1) == ["pair_1", "pair_2", "pair_3"]
    assert read_fastq_names(out2) == ["pair_1", "pair_2", "pair_3"]


def test_pair_preserve_mode_renames_both_mates_consistently(tmp_path):
    r1, r2 = tmp_path / "r1.fq", tmp_path / "r2.fq"
    out1, out2 = tmp_path / "out1.fq.gz", tmp_path / "out2.fq.gz"
    write_fastq(r1, [("same/1", "ACGT"), ("same/1", "AAAA")])
    write_fastq(r2, [("same/2", "TGCA"), ("same/2", "TTTT")])

    fastx_rename_pair([str(r1)], [str(r2)], str(out1), str(out2), mode="preserve")

    assert read_fastq_names(out1) == ["same/1", "same_2/1"]
    assert read_fastq_names(out2) == ["same/2", "same_2/2"]


def test_pair_preserve_avoids_generated_suffix_collisions(tmp_path):
    r1, r2 = tmp_path / "r1.fq", tmp_path / "r2.fq"
    out1, out2 = tmp_path / "out1.fq.gz", tmp_path / "out2.fq.gz"
    names = ["readA", "readA", "readA_2", "readA"]
    write_fastq(r1, [(name + "/1", "ACGT") for name in names])
    write_fastq(r2, [(name + "/2", "TGCA") for name in names])

    fastx_rename_pair([str(r1)], [str(r2)], str(out1), str(out2), mode="preserve")

    assert read_fastq_names(out1) == ["readA/1", "readA_2/1", "readA_2_2/1", "readA_3/1"]
    assert read_fastq_names(out2) == ["readA/2", "readA_2/2", "readA_2_2/2", "readA_3/2"]


def test_pair_mode_rejects_mismatched_file_lists(tmp_path):
    with pytest.raises(ValueError, match="same non-zero number"):
        fastx_rename_pair(["r1a", "r1b"], ["r2a"], str(tmp_path / "o1"), str(tmp_path / "o2"))


def test_gzip_defaults_and_custom_level(tmp_path):
    assert inspect.signature(OpenFqGzip).parameters["compresslevel"].default == 4
    input_fq = tmp_path / "input.fq"
    output_fq = tmp_path / "output.fq.gz"
    write_fastq(input_fq, [("a", "ACGT")])

    fastx_rename(str(input_fq), str(output_fq), mode="index", compresslevel=1)

    assert output_fq.read_bytes()[8] == 4
