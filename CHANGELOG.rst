=========
Changelog
=========

Version 1.3.6 (2026-08-27)
==========================

Command-line interface
----------------------

- Added the top-level ``pybiotk`` command, which lists every console tool
  registered by the installed package together with a short description.
- Added automatic console-entry discovery so the list stays synchronized with
  installed package metadata.
- Expanded command documentation and added tests that require every registered
  tool to have a description.

Performance and utilities
-------------------------

- Reimplemented ``bed2bedgraph`` with an interval-endpoint sweep instead of
  per-base expansion, with support for unsorted, compressed, and multi-file
  BED input.
- Corrected half-open bedGraph coordinates and documented the equivalent
  ``bedtools genomecov`` workflow and chromosome-size generation.
- Improved the bundled FASTA, FASTQ, GTF, chromosome-size, subsequence, and
  duplicate-removal scripts with safer temporary files, input validation,
  gzip support, and portable shell behavior.
- Streamed SAM output in ``merge_subseq.sh`` and fixed the
  ``subseq_analysis`` command-line argument forwarding error.

FASTX and table processing
--------------------------

- Made compact indexed names the default in ``fastx_rename`` while preserving
  matching paired-end identifiers and input order.
- Fixed preserve-mode collisions for repeated read names and strengthened
  paired-end name and record-count validation.
- Fixed ``read_tables`` name filtering, reusable direct API behavior, and
  non-interactive standard-input handling.
- Removed a pandas ``groupby.apply`` deprecation warning from ``merge_row``.

Version 1.3.5 (2026-07-22)
==========================

FASTQ processing
----------------

- Added paired-end mode to ``fastx_rename`` while retaining the existing
  single-end mode.
- Added ``index`` and ``preserve`` rename modes. ``index`` is now the default
  and supports compact base-10, hexadecimal, or base-36 identifiers.
- Preserved read order and assigned matching identifiers to paired-end mates,
  including across multiple input lanes.
- Changed ``fastq_uniq`` to use stable 128-bit BLAKE2b digest keys by default,
  with an ``exact`` comparison mode for collision-free deduplication.
- Corrected duplicate statistics to report duplication rate among eligible
  reads after length filtering.
- Added stricter paired-end validation for unequal record counts and mismatched
  read names, including ``/1`` and ``/2`` mate suffixes.
- Added configurable gzip compression levels to ``fastx_rename`` and
  ``fastq_uniq``. The default compression level is now 4.

Logging
-------

- Added centralized logging configuration through ``get_logger`` and
  ``configure_logging``.
- Added standard ``StreamHandler`` output for library use and Rich-formatted
  stderr output for command-line tools.
- Retained compatibility with ``from pybiotk.utils import logging``.
- Migrated command-line modules to named loggers without repeatedly defining
  logging configuration.

Compatibility and packaging
---------------------------

- Added default arguments to ``read_tables.main`` for direct API use.
- Switched package version resolution to Git tags through ``setuptools_scm``.
- Fixed a Cython interval-tree random-priority edge case and removed the
  associated compiler warning.
- Added tests for FASTQ deduplication, FASTX renaming, paired-end validation,
  gzip configuration, logging behavior, and stream utilities.
