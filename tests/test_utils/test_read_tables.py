from io import StringIO

import pandas as pd
import pytest

from pybiotk.utils.read_tables import main, run


def write_table(path, rows):
    path.write_text("name\tvalue\n" + "".join(f"{name}\t{value}\n" for name, value in rows))


def read_output(path):
    return pd.read_table(path, dtype=str)


def test_name_path_filters_every_input_table(tmp_path):
    first = tmp_path / "first.tsv"
    second = tmp_path / "second.tsv"
    names = tmp_path / "names.txt"
    output = tmp_path / "output.tsv"
    write_table(first, [("keep", "1"), ("drop", "2")])
    write_table(second, [("drop", "3"), ("keep", "4")])
    names.write_text("keep\n")

    main([str(first), str(second)], str(output), str(names))

    assert read_output(output).to_dict("records") == [
        {"name": "keep", "value": "1"},
        {"name": "keep", "value": "4"},
    ]


def test_name_stream_is_loaded_once_for_multiple_tables(tmp_path):
    first = tmp_path / "first.tsv"
    second = tmp_path / "second.tsv"
    output = tmp_path / "output.tsv"
    write_table(first, [("keep", "1"), ("drop", "2")])
    write_table(second, [("drop", "3"), ("keep", "4")])

    main([str(first), str(second)], str(output), StringIO("keep\n"))

    assert read_output(output)["value"].tolist() == ["1", "4"]


def test_contains_treats_names_as_literal_substrings_and_handles_missing_values(tmp_path):
    input_table = tmp_path / "input.tsv"
    names = tmp_path / "names.txt"
    output = tmp_path / "output.tsv"
    write_table(input_table, [("prefix-a.b-suffix", "1"), ("axb", "2"), ("", "3")])
    names.write_text("a.b\n")

    main([str(input_table)], str(output), str(names), contains=True)

    assert read_output(output)["value"].tolist() == ["1"]


@pytest.mark.parametrize(
    ("exclude", "expected"),
    [(False, []), (True, ["1", "2"])],
)
def test_empty_name_list_has_set_filter_semantics(tmp_path, exclude, expected):
    input_table = tmp_path / "input.tsv"
    names = tmp_path / "names.txt"
    output = tmp_path / "output.tsv"
    write_table(input_table, [("one", "1"), ("two", "2")])
    names.write_text("")

    main([str(input_table)], str(output), str(names), exclude=exclude)

    assert read_output(output)["value"].tolist() == expected


def test_delimiter_is_literal_and_missing_values_do_not_fail(tmp_path):
    input_table = tmp_path / "input.tsv"
    names = tmp_path / "names.txt"
    output = tmp_path / "output.tsv"
    write_table(input_table, [("one.two", "1"), ("one-two", "2"), ("", "3")])
    names.write_text("two\n")

    main([str(input_table)], str(output), str(names), delimiter=".")

    assert read_output(output)["value"].tolist() == ["1"]


def test_invalid_column_has_clear_error(tmp_path):
    input_table = tmp_path / "input.tsv"
    names = tmp_path / "names.txt"
    output = tmp_path / "output.tsv"
    write_table(input_table, [("one", "1")])
    names.write_text("one\n")

    with pytest.raises(ValueError, match="column index 2 is out of range"):
        main([str(input_table)], str(output), str(names), column=2)


def test_non_tty_stdin_does_not_enable_name_filtering(tmp_path, monkeypatch):
    input_table = tmp_path / "input.tsv"
    output = tmp_path / "output.tsv"
    write_table(input_table, [("one", "1"), ("two", "2")])
    monkeypatch.setattr("sys.stdin", StringIO("one\n"))
    monkeypatch.setattr("sys.argv", ["read_tables", str(input_table), "-o", str(output)])

    run()

    assert read_output(output)["value"].tolist() == ["1", "2"]


def test_explicit_stdin_name_filter_is_preserved(tmp_path, monkeypatch):
    input_table = tmp_path / "input.tsv"
    output = tmp_path / "output.tsv"
    write_table(input_table, [("one", "1"), ("two", "2")])
    monkeypatch.setattr("sys.stdin", StringIO("one\n"))
    monkeypatch.setattr(
        "sys.argv", ["read_tables", str(input_table), "-o", str(output), "-n", "-"]
    )

    run()

    assert read_output(output)["value"].tolist() == ["1"]
