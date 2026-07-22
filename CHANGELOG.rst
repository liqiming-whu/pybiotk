=========
Changelog
=========

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
