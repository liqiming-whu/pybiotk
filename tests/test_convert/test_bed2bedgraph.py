import gzip

import pytest

from pybiotk.convert.bed2bedgraph import main


def write_bed(path, lines):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt") as stream:
        stream.writelines(lines)


def test_bed2bedgraph_uses_half_open_coverage_segments(tmp_path, capsys):
    bed = tmp_path / "reads.bed"
    write_bed(bed, ["chr1\t2\t5\n", "chr1\t0\t3\n", "chr1\t5\t6\n"])

    main([str(bed)], False)

    assert capsys.readouterr().out == (
        "chr1\t0\t2\t1\n"
        "chr1\t2\t3\t2\n"
        "chr1\t3\t6\t1\n"
    )


def test_bed2bedgraph_combines_files_and_sorts_chromosomes(tmp_path, capsys):
    bed1 = tmp_path / "lane1.bed"
    bed2 = tmp_path / "lane2.bed.gz"
    write_bed(bed1, ["track name=test\n", "chr2\t0\t2\n", "chr1\t1\t2\n"])
    write_bed(bed2, ["# comment\n", "\n", "chr1\t0\t1\n"])

    main([str(bed1), str(bed2)], False)

    assert capsys.readouterr().out == "chr1\t0\t2\t1\nchr2\t0\t2\t1\n"


def test_bed2bedgraph_header_uses_clean_filename(tmp_path, capsys):
    bed = tmp_path / "sample.bed.gz"
    write_bed(bed, ["chr1\t0\t1\n"])

    main([str(bed)], True)

    output = capsys.readouterr().out.splitlines()
    assert 'name="sample"' in output[0]
    assert output[1] == "chr1\t0\t1\t1"


@pytest.mark.parametrize("line", ["chr1\t0\n", "chr1\tx\t2\n", "chr1\t2\t2\n"])
def test_bed2bedgraph_rejects_invalid_records(tmp_path, line):
    bed = tmp_path / "invalid.bed"
    write_bed(bed, [line])

    with pytest.raises(ValueError):
        main([str(bed)], False)
