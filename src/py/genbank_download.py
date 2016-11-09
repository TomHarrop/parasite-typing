#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from Bio import Entrez
from Bio import SeqIO


###############
# FUNCTIONS ###
###############


# get entrez GIs for COI genes for species
def get_coi_gis(species_name):
    term = (species_name + '[ORGN] AND '
            '(COI OR cytochrome oxidase subunit I)')
    handle = Entrez.esearch(db='nuccore', term=term)
    record = Entrez.read(handle)
    gi_list = record['IdList']
    return(gi_list)


# get entrez GIs for ITS1 genes for species
def get_its_gis(species_name):
    term = (species_name + '[ORGN] AND '
            '(ITS1 OR internal transcribed spacer)')
    handle = Entrez.esearch(db='nuccore', term=term)
    record = Entrez.read(handle)
    gi_list = record['IdList']
    return(gi_list)


# flatten list cf. `unlist()` in R
def flatten_list(l):
    for x in l:
        if hasattr(x, '__iter__') and not isinstance(x, str):
            for y in flatten_list(x):
                yield y
        else:
            yield x

########
# CODE #
########

def main():

    # require email address from cli
    parser = argparse.ArgumentParser()
    parser.add_argument('-e',
                        help='email address',
                        metavar='email',
                        type=str,
                        required=True)
    args = parser.parse_args()

    # identify myself to Entrez
    Entrez.email = args.e

    # which fields can I search?
    handle = Entrez.einfo()
    record = Entrez.read(handle)
    handle = Entrez.einfo(db='nuccore')
    record = Entrez.read(handle)
    record['DbInfo']['FieldList']

    # Search for COI genes in weevil and both wasps
    species_names = ['microctonus hyperodae', 'microctonus aethiopoides',
                     'listronotus bonariensis']
    gi_list_nested = list(get_coi_gis(x) for x in species_names)
    gi_list_all = list(flatten_list(gi_list_nested))

    # download records
    handle = Entrez.efetch(
        db='nuccore',
        id=gi_list_all,
        rettype='gb',
        retmode='text')
    coi_records = SeqIO.parse(handle, 'gb')

    # output FASTA
    SeqIO.write(coi_records, 'data/coi.fa', 'fasta')

    # get ITS gis
    its_gis_nested = list(get_its_gis(x) for x in species_names)
    its_gis_all = list(flatten_list(its_gis_nested))

    # download ITS records
    handle = Entrez.efetch(
        db='nuccore',
        id=its_gis_all,
        rettype='gb',
        retmode='text')
    its_records = SeqIO.parse(handle, 'gb')

    # write fasta
    SeqIO.write(its_records, 'data/its.fa', 'fasta')

if __name__ == '__main__':
    main()
