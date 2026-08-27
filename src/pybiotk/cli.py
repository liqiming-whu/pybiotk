"""List command-line tools installed with pybiotk."""
import argparse
import sys
from importlib.metadata import PackageNotFoundError, distribution
from typing import Dict, List, Sequence, Tuple


TOOL_DESCRIPTIONS: Dict[str, str] = {
    "bam2fastx": "Convert BAM alignments to FASTA or FASTQ.",
    "bam_random": "Randomly sample alignments from a BAM file.",
    "bampe_order_by_name": "Write paired-end BAM alignments in query-name order.",
    "bed2bedgraph": "Convert BED intervals to bedGraph coverage.",
    "bigwigfetcher": "Fetch signal values from bigWig files.",
    "count_normalize": "Normalize count-table columns.",
    "ercc_parser": "Convert an ERCC reference table to FASTA and GTF.",
    "fa2fastq": "Convert FASTA records to FASTQ.",
    "fasta_filter": "Filter FASTA records by name or sequence length.",
    "fastq_join": "Join overlapping paired-end FASTQ reads.",
    "fastq_uniq": "Remove duplicate FASTQ reads or read pairs.",
    "fastx_rename": "Rename single-end or paired-end FASTA/FASTQ records.",
    "fq2fasta": "Convert FASTQ records to FASTA.",
    "genomefetcher": "Fetch genomic sequences for regions or annotations.",
    "gtf2bed": "Convert GTF annotations to BED.",
    "gtf_filter": "Filter GTF records by feature or attribute.",
    "infer_experiment": "Infer RNA-seq layout and strandedness from alignments.",
    "merge_row": "Merge duplicate table rows using selected columns.",
    "merge_subseq": "Merge FASTA sequences contained within other sequences.",
    "merge_transcript": "Merge transcript annotations by genomic overlap.",
    "metaplot": "Create aggregate signal profiles around genomic features.",
    "pyanno": "Annotate BED or BAM intervals with GTF features.",
    "read_tables": "Join, select, and filter tabular files.",
    "reference_count": "Count BAM reads or fragments by reference sequence.",
    "reverse_fastx": "Reverse or reverse-complement FASTA/FASTQ sequences.",
    "rmats_filter": "Filter and rank rMATS alternative-splicing results.",
    "rna_fragment_size": "Calculate RNA fragment-length distributions from BAM.",
    "seq_random": "Generate random genomic sequences.",
    "subseq_analysis": "Analyze sequence-containment relationships.",
    "summary_log": "Extract metrics from bioinformatics log files.",
}


def registered_tools() -> List[str]:
    """Return console scripts registered by the installed pybiotk package."""
    try:
        entry_points = distribution("pybiotk").entry_points
    except PackageNotFoundError:
        return sorted(TOOL_DESCRIPTIONS)
    return sorted(
        entry_point.name
        for entry_point in entry_points
        if entry_point.group == "console_scripts" and entry_point.name != "pybiotk"
    )


def tool_rows(tool_names: Sequence[str]) -> List[Tuple[str, str]]:
    """Return tool names paired with short descriptions."""
    return [
        (name, TOOL_DESCRIPTIONS.get(name, "Run this command with --help for usage details."))
        for name in tool_names
    ]


def format_tool_list(rows: Sequence[Tuple[str, str]]) -> str:
    """Format tool descriptions as a terminal-friendly plain-text table."""
    if not rows:
        return "No pybiotk command-line tools are registered."
    width = max(len(name) for name, _ in rows)
    lines = ["Available pybiotk tools:", ""]
    lines.extend(f"  {name:<{width}}  {description}" for name, description in rows)
    lines.extend(["", "Run '<tool> --help' for detailed usage."])
    return "\n".join(lines)


def run(argv: Sequence[str] = None) -> None:
    parser = argparse.ArgumentParser(
        prog="pybiotk",
        description="List command-line tools installed with pybiotk.",
    )
    parser.parse_args(argv)
    print(format_tool_list(tool_rows(registered_tools())))


if __name__ == "__main__":
    run(sys.argv[1:])
