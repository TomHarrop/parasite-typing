#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import SeqIO
import re
import os
import tompytools


def rename_records(records):
    for rec in records:
        re_spec = re.compile(r'^(\w)\S*\s(\w)')
        re_match = re_spec.search(rec.description)
        spec_sn = re_match.group(1) + re_match.group(2)
        new_id = spec_sn + "_" + rec.id
        rec.id = new_id
    return(records)


def main():

    # scan data directory for GenBank files
    genbank_files = {os.path.splitext(x.name)[0]: x.path
                     for x in os.scandir('data') if
                     x.name.endswith('.gb') and x.is_file}

    # convert each file to fasta
    for file_name in genbank_files:
        tompytools.generate_message('Reading GenBank from ' +
                                    genbank_files[file_name])
        # parse GenBank
        with open(genbank_files[file_name], 'r') as handle:
            gb_records = list(SeqIO.parse(handle, 'gb'))

        # rename and convert to FASTA
        tompytools.generate_message('Renaming')
        renamed_records = rename_records(gb_records)

        # write output
        output_file = 'output/fa/' + file_name + '.fa'
        tompytools.generate_message('Writing output to ' + output_file)
        SeqIO.write(renamed_records, output_file, "fasta")

if __name__ == '__main__':
    main()
