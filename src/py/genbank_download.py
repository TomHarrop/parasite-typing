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
def get_gis(species_name, search_term):
    term = (species_name + '[ORGN] AND ' + search_term)
    print('Search term: ' + term)
    handle = Entrez.esearch(db='nuccore', term=term)
    record = Entrez.read(handle)
    gi_list = record['IdList']
    return(gi_list)


# get GenBank records for GIs
def get_gb_records(gi_list):
    handle = Entrez.efetch(
        db='nuccore',
        id=gi_list,
        rettype='gb',
        retmode='text')
    return(SeqIO.parse(handle, 'gb'))


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
    species_names = {'mh': 'microctonus hyperodae',
                     'ma': 'microctonus aethiopoides',
                     'lb': 'listronotus bonariensis'}

    # genes to search for
    search_terms = {'coi': '(COI OR cytochrome oxidase subunit I)',
                    'its1': '(ITS1 OR internal transcribed spacer)',
                    '16s': '16S',
                    '28s': '28S'}

    # start seearch
    for term_key in search_terms:
        tompytools.generate_message('Getting GIs for ' + term_key)

        # search by species for each search_term
        for spec in species_names:
            tompytools.generate_message('Running search for ' + spec)
            gi_results = get_gis(species_name=species_names[spec],
                                 search_term=search_terms[term_key])

            # download results
            tompytools.generate_message('Downloading GenBank records')
            gb_records = get_gb_records(gi_results)

            # write gb data to disk
            file_name = 'data/' + spec + '_' + term_key + '.gb'
            tompytools.generate_message('Writing GenBank to ' + file_name)
            SeqIO.write(gb_records, file_name, 'gb')

    tompytools.generate_message('Done')

if __name__ == '__main__':
    main()
