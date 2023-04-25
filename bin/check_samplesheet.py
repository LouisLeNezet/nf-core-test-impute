#!/usr/bin/env python


"""Provide a command line tool to validate and transform tabular samplesheets."""


import argparse
import csv
import logging
import sys
from collections import Counter
from pathlib import Path

logger = logging.getLogger()


class RowChecker:
    """
    Define a service that can validate and transform each given row.

    Attributes:
        modified (list): A list of dicts, where each dict corresponds to a previously
            validated and transformed row. The order of rows is maintained.

    """

    VALID_BAM_FORMATS = (
        ".bam",
        ".cram"
    )
    VALID_BAI_FORMATS = (
        ".bai",
        ".crai"
    )

    def __init__(
        self,
        sample_col="sample_id",
        bam_col="bam",
        bai_col="bam_index",
        ref_col="ref_id",
        **kwargs,
    ):
        """
        Initialize the row checker with the expected column names.

        Args:
            sample_col (str): The name of the column that contains the sample name
                (default "sample_id").
            bam_col (str): The name of the column that contains the bam file path (default "bam").
            bai_col (str): The name of the column that contains the bai file path (default "bai").
            ref_col (str): The name of the column that contains the reference genome id (default "ref_id").

        """
        super().__init__(**kwargs)
        self._sample_col = sample_col
        self._bam_col = bam_col
        self._bai_col = bai_col
        self._ref_col = ref_col
        self._seen = set()
        self.modified = []

    def validate_and_transform(self, row):
        """
        Perform all validations on the given row and insert the read pairing status.

        Args:
            row (dict): A mapping from column headers (keys) to elements of that row
                (values).

        """
        self._validate_sample(row)
        self._validate_bam(row)
        self._validate_bai(row)
        self._validate_ref(row)
        self._seen.add((row[self._sample_col], row[self._bam_col]))
        self.modified.append(row)

    def _validate_sample(self, row):
        """Assert that the sample name exists and convert spaces to underscores."""
        if len(row[self._sample_col]) <= 0:
            raise AssertionError("Sample input is required.")
        # Sanitize samples slightly.
        row[self._sample_col] = row[self._sample_col].replace(" ", "_")

    def _validate_bam(self, row):
        """Assert that the BAM entry is non-empty and has the right format."""
        if len(row[self._bam_col]) <= 0:
            raise AssertionError("The BAM file is required.")
        self._validate_bam_format(row[self._bam_col])

    def _validate_bai(self, row):
        """Assert that the BAI entry is non-empty and has the right format."""
        if len(row[self._bai_col]) <= 0:
            raise AssertionError("The BAI file is required.")
        self._validate_bai_format(row[self._bai_col])

    def _validate_ref(self, row):
        """Assert that the ref name exists and convert spaces to underscores."""
        if len(row[self._ref_col]) <= 0:
            raise AssertionError("Ref input is required.")
        # Sanitize ref slightly.
        row[self._ref_col] = row[self._ref_col].replace(" ", "_")

    def _validate_bam_format(self, filename):
        """Assert that a given filename has one of the expected bam extensions."""
        if not any(filename.endswith(extension) for extension in self.VALID_BAM_FORMATS):
            raise AssertionError(
                f"The BAM file has an unrecognized extension: {filename}\n"
                f"It should be one of: {', '.join(self.VALID_BAM_FORMATS)}"
            )

    def _validate_bai_format(self, filename):
        """Assert that a given filename has one of the expected bam extensions."""
        if not any(filename.endswith(extension) for extension in self.VALID_BAI_FORMATS):
            raise AssertionError(
                f"The BAI file has an unrecognized extension: {filename}\n"
                f"It should be one of: {', '.join(self.VALID_BAI_FORMATS)}"
            )
            

    def validate_unique_samples(self):
        """
        Assert that the combination of sample name and BAM filename is unique.

        In addition to the validation, also rename all samples to have a suffix of _T{n}, where n is the
        number of times the same sample exist, but with different BAM files, e.g., multiple runs per experiment.

        """
        if len(self._seen) != len(self.modified):
            raise AssertionError("The pair of sample name and BAM must be unique.")
        seen = Counter()
        for row in self.modified:
            sample = row[self._sample_col]
            seen[sample] += 1
            row[self._sample_col] = f"{sample}_T{seen[sample]}"


def read_head(handle, num_lines=10):
    """Read the specified number of lines from the current position in the file."""
    lines = []
    for idx, line in enumerate(handle):
        if idx == num_lines:
            break
        lines.append(line)
    return "".join(lines)


def sniff_format(handle):
    """
    Detect the tabular format.

    Args:
        handle (text file): A handle to a `text file`_ object. The read position is
        expected to be at the beginning (index 0).

    Returns:
        csv.Dialect: The detected tabular format.

    .. _text file:
        https://docs.python.org/3/glossary.html#term-text-file

    """
    peek = read_head(handle)
    handle.seek(0)
    sniffer = csv.Sniffer()
    if not sniffer.has_header(peek):
        logger.critical("The given sample sheet does not appear to contain a header.")
        sys.exit(1)
    dialect = sniffer.sniff(peek)
    return dialect


def check_samplesheet(file_in, file_out):
    """
    Check that the tabular samplesheet has the structure expected by nf-core pipelines.

    Validate the general shape of the table, expected columns, and each row. Also add
    an additional column which records the BAM file and the reference genome ID.

    Args:
        file_in (pathlib.Path): The given tabular samplesheet. The format can be either
            CSV, TSV, or any other format automatically recognized by ``csv.Sniffer``.
        file_out (pathlib.Path): Where the validated and transformed samplesheet should
            be created; always in CSV format.

    Example:
        This function checks that the samplesheet follows the following structure,
        see also the `viral recon samplesheet`_::

            sample_id,bam,bam_index,ref_id
            SAMPLE_PE,SAMPLE_PE_RUN1_BAM,SAMPLE_PE_RUN1_BAM_BAI,REF_1
            SAMPLE_PE,SAMPLE_PE_RUN2_BAM,SAMPLE_PE_RUN2_BAM_BAI,REF_1

    .. _viral recon samplesheet:
        https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv

    """
    required_columns = {"sample_id", "bam", "bam_index", "ref_id"}
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_in.open(newline="") as in_handle:
        reader = csv.DictReader(in_handle, dialect=sniff_format(in_handle))
        # Validate the existence of the expected header columns.
        if not required_columns.issubset(reader.fieldnames):
            req_cols = ", ".join(required_columns)
            logger.critical(f"The sample sheet **must** contain these column headers: {req_cols}.")
            sys.exit(1)
        # Validate each row.
        checker = RowChecker()
        for i, row in enumerate(reader):
            try:
                checker.validate_and_transform(row)
            except AssertionError as error:
                logger.critical(f"{str(error)} On line {i + 2}.")
                sys.exit(1)
        checker.validate_unique_samples()
    header = list(reader.fieldnames)
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_out.open(mode="w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, header, delimiter=",")
        writer.writeheader()
        for row in checker.modified:
            writer.writerow(row)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv samplesheet.valid.csv",
    )
    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="Tabular input samplesheet in CSV or TSV format.",
    )
    parser.add_argument(
        "file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Transformed output samplesheet in CSV format.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    check_samplesheet(args.file_in, args.file_out)


if __name__ == "__main__":
    sys.exit(main())
