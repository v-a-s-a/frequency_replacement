#!/usr/bin/env python
"""
A frequency replacement program for DANER files.
"""

import gzip
import datetime
import argparse
import sys
import numpy as np


def smart_open(f):
    if f.endswith(".gz"):
        fconn = gzip.open(f, 'r')
    else:
        fconn = open(f, 'r')
    return fconn


class Variant:

    def __init__(self, ID, chromosome, position, allele_1, allele_2,
                 allele_1_frequency, allele_2_frequency=None):
        self.ID = ID
        self.chromosome = chromosome
        self.position = position
        self.allele_1 = allele_1
        self.allele_2 = allele_2
        self.allele_1_frequency = allele_1_frequency

        if allele_2_frequency:
            self.allele_2_frequency = allele_2_frequency
        else:
            self.allele_2_frequency = str(1.0 - float(allele_1_frequency))

        self.alleles = set((allele_1, allele_2))

    @classmethod
    def from_daner_row(cls, row, header_index, A1_FRQ_field="FRQ_A",
                       A2_FRQ_field=None):
        A1_FRQ = next(filter(lambda x: x.startswith(A1_FRQ_field),
                             header_index.keys()))
        variant = cls(row[header_index["SNP"]],
                      row[header_index["CHR"]],
                      row[header_index["BP"]],
                      row[header_index["A1"]],
                      row[header_index["A2"]],
                      row[header_index[A1_FRQ]])
        if A2_FRQ_field:
            A2_FRQ = next(filter(lambda x: x.startswith(A2_FRQ_field),
                                 header_index.keys()))
            variant = cls(row[header_index["SNP"]],
                          row[header_index["CHR"]],
                          row[header_index["BP"]],
                          row[header_index["A1"]],
                          row[header_index["A2"]],
                          row[header_index[A1_FRQ]],
                          row[header_index[A2_FRQ]])
        return variant

    def __eq__(self, other):
        return self.ID == other.ID and self.chromosome == other.chromosome and self.position == other.position and self.alleles == other.alleles

    def __hash__(self):
        return hash(str(self))

    def __repr__(self):
        return '\t'.join(
            [self.ID, self.chromosome, self.position, str(self.alleles)])


class Daner:

    def __init__(self, daner_file, reference_files):
        self.daner_file = daner_file

        with smart_open(self.daner_file) as fconn:
            self.header = next(fconn).decode("utf-8").strip()
            self.header_index = {field: index for index,
                                 field in enumerate(self.header.split())}

        ncase_header = next(
            filter(lambda x: x.startswith("FRQ_A_"), self.header.split()))
        ncontrol_header = next(
            filter(lambda x: x.startswith("FRQ_U_"), self.header.split()))
        self.n_case = ncase_header.lstrip("FRQ_A_")
        self.n_control = ncontrol_header.lstrip("FRQ_U_")

        self.variants_filtered = list()
        self.variants_not_in_reference = list()
        self.total_variants = 0
        self.frequency_differences = set()

        self.keep_until_column = 12

        self.reference = Reference(reference_files)

    def print_header(self):
        return '\t'.join(
            self.header.split()[:self.keep_until_column])

    def parse_row(self, row):
        row = row.split()
        self.total_variants += 1
        variant = Variant.from_daner_row(row, self.header_index)
        # apply row filters
        if row[self.header_index["Direction"]].count('?') >= 2:
            self.variants_filtered.append('\t'.join(row))
            row = None
        else:
            # replace allele frequency
            frequency = self.reference.get_frequency(variant)
            if frequency:
                self.frequency_differences.add(
                    float(self.header_index["FRQ_U_" + self.n_control]) - float(frequency))
                row[self.header_index[
                    "FRQ_A_{}".format(self.n_case)]] = frequency
                row[self.header_index["FRQ_U_{}".format(
                    self.n_control)]] = frequency
                row = '\t'.join(row[:self.keep_until_column])
            else:
                self.variants_not_in_reference.append('\t'.join(row))
                row = None
        return row

    def replacement_row_generator(self):
        with smart_open(self.daner_file) as fconn:
            next(fconn)
            for row in fconn:
                row = row.decode("utf-8")
                new_row = self.parse_row(row)
                if new_row:
                    yield new_row
                else:
                    continue


class Reference:

    def __init__(self, reference_file_list):
        self.reference_file_list = reference_file_list
        self.variants = dict()

        self.header_index = {"SNP": 0, "A1": 2, "A1_FRQ": 3, "A2": 4,
                             "A2_FRQ": 5, "CHR": 7, "BP": 8}

        for reference_file in self.reference_file_list:
            with smart_open(reference_file) as fconn:
                for row in fconn:
                    row = row.split()
                    variant = Variant.from_daner_row(row,
                                                     self.header_index,
                                                     A1_FRQ_field="A1_FRQ",
                                                     A2_FRQ_field="A2_FRQ")
                    self.variants[variant] = variant

    def get_frequency(self, variant):
        reference_variant = self.variants.get(variant)
        if reference_variant:
            if variant.allele_1 is reference_variant.allele_1:
                reference_frequency = reference_variant.allele_1_frequency
            elif variant.allele_1 is reference_variant.allele_2:
                reference_frequency = reference_variant.allele_2_frequency
        else:
            reference_frequency = None
        return reference_frequency


def __main__():

    parser = argparse.ArgumentParser(
        description="Replacing frequencies from a DANER file given a reference.")
    parser.add_argument("--out", action="store", dest="out_daner",
                        help="Output location of the resulting gzipped file.")
    parser.add_argument("--readme-file", action="store", dest="out_readme",
                        help="Output location of the accompanying README file.")
    parser.add_argument("--references", action="store", nargs='+',
                        dest="references", help="List of reference file locations used to replace frequencies.")
    parser.add_argument("--daner", action="store", dest="daner",
                        help="Location of gzipped DANER file to replace frequencies in.")
    args = parser.parse_args()

    # Generate a README for this replacement
    readme_template = '''
# {{DATASET NAME}} Public Release

## Pre-release Filters
In addition to study level and meta-analysis quality control, we apply a final
variant filter. This currently consists of:
  - Remove variants that are missing in two or more studies.

## Frequency Replacement
Allele frequency estimates from the study are replaced with estimates from
{{AUTOSOMAL REFERENCE NAME}} and {{CHR23 REFERENCE NAME}} to reduce the ability
 for a third party to reidentify study participants.
The frequency replacement proceeds as follows:
  1. Load reference European frequencies from {{AUTOSOMAL REFERENCE NAME}} for
        autosomal variants.
  2. Load reference European frequencies from {{CHR23 REFERENCE NAME}} release
        for chr23 variants.
  3. Merge reference and study variants on variant ID, chromosome and physical
        position.
  4. Invert reference allele frequencies if the reference allele codes are
        flipped.
  5. Leave variant allele frequencies that are not found in the reference (a
        very small proportion of variants).
  6. Set both case and control frequency estimates to the same value.

## Misc.
Command used to generate files:
{command_line}


Last updated: {todays_date}

    '''
    readme_text = readme_template.format(
        todays_date=datetime.date.today().strftime("%B %d, %Y"),
        command_line=' '.join(sys.argv))

    # write readme to file
    with open(args.out_readme, 'w') as readme_conn:
        readme_conn.write(readme_text)

    # write DANER with replaced frequencies to file
    pgc_data = Daner(
        daner_file=args.daner,
        reference_files=args.references)
    with gzip.open(args.out_daner, 'wb') as daner_conn:
        daner_conn.write(bytes(pgc_data.print_header() + '\n', "utf-8"))
        for row in pgc_data.replacement_row_generator():
            daner_conn.write(bytes(row + '\n', "utf-8"))

    # write all missing variants to file
    with gzip.open(args.out_daner.rstrip(".gz") + ".missing.gz", 'wb') as missing_conn:
        missing_conn.write(bytes(pgc_data.print_header() + '\n', "utf-8"))
        for variant in pgc_data.variants_not_in_reference:
            missing_conn.write(bytes(variant + '\n', "utf-8"))

    # write all filtered variants to file
    with gzip.open(args.out_daner.rstrip(".gz") + ".filtered.gz", 'wb') as filtered_conn:
        filtered_conn.write(bytes(pgc_data.header + '\n', "utf-8"))
        for variant in pgc_data.variants_filtered:
            filtered_conn.write(bytes(variant + '\n', "utf-8"))

    summary_message = """
Summary:
    {var_in_pgc} variants in orginal DANER file.
    {var_filtered} variants filtered from original DANER.
    {var_in_references} variants in all references.
    {var_not_found} variants not found in references.

    {extreme_freq_diffs} variants with >10% change in frequency.
    """

    print(summary_message.format(
        var_in_pgc=pgc_data.total_variants,
        var_in_references=len(pgc_data.reference.variants),
        var_not_found=len(pgc_data.variants_not_in_reference),
        var_filtered=len(pgc_data.variants_filtered),
        extreme_freq_diffs=sum(
            map(lambda x: np.abs(x) >= 0.1, pgc_data.frequency_differences))))


if __name__ == "__main__":
    __main__()
