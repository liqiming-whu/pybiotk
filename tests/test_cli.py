from pathlib import Path
import tomllib

from pybiotk import cli


def configured_tools():
    with (Path(__file__).parents[1] / "pyproject.toml").open("rb") as handle:
        project = tomllib.load(handle)["project"]
    return set(project["scripts"])


def test_descriptions_cover_all_registered_tools():
    assert configured_tools() == {"pybiotk", *cli.TOOL_DESCRIPTIONS}


def test_tool_rows_use_fallback_description():
    assert cli.tool_rows(["future_tool"]) == [
        ("future_tool", "Run this command with --help for usage details.")
    ]


def test_run_prints_tools(monkeypatch, capsys):
    monkeypatch.setattr(cli, "registered_tools", lambda: ["gtf2bed", "pyanno"])

    cli.run([])

    output = capsys.readouterr().out
    assert "Available pybiotk tools:" in output
    assert "gtf2bed" in output
    assert "pyanno" in output
    assert "<tool> --help" in output
