#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from Bio import Entrez
from Bio import SeqIO
import tompytools

###############
# FUNCTIONS ###
###############


# get entrez GIs for COI genes for species
def get_coi_gis(species_name):
    term = (species_name + '[ORGN] AND '
            '(COI OR cytochrome oxidase subunit I)')
    print('Search term: ' + term)
    handle = Entrez.esearch(db='nuccore', term=term)
    record = Entrez.read(handle)
    gi_list = record['IdList']
    return(gi_list)


# get entrez GIs for ITS1 genes for species
def get_its_gis(species_name):
    term = (species_name + '[ORGN] AND '
            '(ITS1 OR internal transcribed spacer)')
    print('Search term: ' + term)
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

    # Search in weevil and both wasps
    species_names = ['microctonus hyperodae', 'microctonus aethiopoides',
                     'listronotus bonariensis']

    # get coi gene ids
    tompytools.generate_message('Getting COI GIs')
    gi_list_nested = list(get_coi_gis(x) for x in species_names)
    gi_list_all = list(flatten_list(gi_list_nested))

    # download records
    tompytools.generate_message('Downloading records')
    handle = Entrez.efetch(
        db='nuccore',
        id=gi_list_all,
        rettype='gb',
        retmode='text')
    coi_records = SeqIO.parse(handle, 'gb')

    # output FASTA
    tompytools.generate_message('Writing FASTA')
    SeqIO.write(coi_records, 'data/coi.gb', 'gb')

    # get ITS gis
    tompytools.generate_message('Getting ITS1 GIs')
    its_gis_nested = list(get_its_gis(x) for x in species_names)
    its_gis_all = list(flatten_list(its_gis_nested))

    # download ITS records
    tompytools.generate_message('Downloading records')
    handle = Entrez.efetch(
        db='nuccore',
        id=its_gis_all,
        rettype='gb',
        retmode='text')
    its_records = SeqIO.parse(handle, 'gb')

    # write fasta
    tompytools.generate_message('Writing FASTA')
    SeqIO.write(its_records, 'data/its.gb', 'gb')

    tompytools.generate_message('Done')

if __name__ == '__main__':
    main()
