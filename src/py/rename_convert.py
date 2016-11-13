#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import SeqIO
import re


def rename_records(records):
    for rec in records:
        re_spec = re.compile(r'^(\w)\S*\s(\w)')
        re_match = re_spec.search(rec.description)
        spec_sn = re_match.group(1) + re_match.group(2)
        new_id = spec_sn + "_" + rec.id
        rec.id = new_id
    return(records)


def main():

    with open('data/coi.gb', 'r') as handle:
            coi_records = list(SeqIO.parse(handle, 'gb'))
    coi_records_renamed = rename_records(coi_records)

    with open('data/its.gb', 'r') as handle:
            its_records = list(SeqIO.parse(handle, 'gb'))
    its_records_renamed = rename_records(its_records)

    SeqIO.write(coi_records_renamed, "output/coi.fa", "fasta")
    SeqIO.write(its_records_renamed, "output/its.fa", "fasta")

if __name__ == '__main__':
    main()
