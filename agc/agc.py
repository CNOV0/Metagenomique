#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
import textwrap
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication  (default 100)")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()

def read_fasta(amplicon_file, minseqlen):
    """
    Read a zipped fasta file
    Parameters
    ----------
    amplicon_file:
    string
        Name of the fasta file
    minseqlen:
    int
        Length of the minimum sequence length
    Returns
    -------
    Generator containing the selected sequences
    """
    isfile(amplicon_file)
    seq  = ""
    with gzip.open(amplicon_file, "rt") as file_in:
        for line in file_in:
            if line.startswith('>'):
                if len(seq) >= minseqlen:
                     yield seq
                seq  = ""
            else:
                seq+=line.strip()
        if len(seq) >= minseqlen:
            yield seq


def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    """
    Sort out the duplicated sequences
    Parameters
    ----------
    amplicon_file:
    string
        Name of the fasta file
    minseqlen:
    int
        Length of the minimum sequence length
    mincount:
    int
        Minimum number of count
    Returns
    -------
    Generator containing the selected sequences and their occurrences
    """
    seq_list = list(read_fasta(amplicon_file, minseqlen))
    seq_uniq = set(seq)
    res = []
    for seq in uniq_seq:
        count = seq.count(seq)
        if count >= mincount:
            res.append([seq, occur])
    res = sorted(res, key=lambda x : x[1], reverse=True)
    for i in res:
        yield i


def get_identity(alignment_list):
    """
    Calculate the percentage of identity between two sequences
    Parameters
    ----------
    amplicon_list:
    list
        Two sequences aligned
    Returns
    -------
    Percentage of identity bewteen the sequences
    """
    nucl_id = 0
    for i in range(len(alignment_list[0])):
        if alignment_list[0][i] == alignment_list[1][i]:
            nucl_id+=1
    return nucl_id / len(alignment_list[0])*100


def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """
    Create a list of OTU with greedy algorithm
    Parameters
    ----------
    amplicon_file:
    string
        Name of the fasta file
    minseqlen:
    int
        Length of the minimum sequence length
    mincount:
    int
        Minimum number of count
    chunk_size:
    int
        Size of the chunk
    kmer_size:
    int
        Size of the wanted kmer
    Returns
    -------
    List of OTU with their occurrence
    """
    seq_list = list(dereplication_fulllength(amplicon_file, minseqlen, mincount))
    OTU_list = [seq_list[0]]
    for seq_1 in seq_list[1:]:
        for seq_2 in OTU_list:
            alignment = nw.global_align(seq_1[0], seq_2[0], gap_open=-1, 
                                             gap_extend=-1, 
                        matrix=os.path.abspath(os.path.join(os.path.dirname(__file__),"MATCH")))
        if get_identity(alignment) < 97:
            OTU_list.append(seq_1)
    return OTU_list

def write_OTU(OTU_list, output_file):
    """
    Calculate the percentage of identity between two sequences
    Parameters
    ----------
    amplicon_list:
    list
        Two sequences aligned
    Returns
    -------
    Percentage of identity bewteen the sequences
    """
    with open(output_file, "w") as file:
        for i, (sequence, count) in enumerate(OTU_list):
            file.write(f">OTU_{i+1} occurrence:{count}\n{textwrap.fill(sequence, width = 80)}\n")

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    
    # Votre programme ici
    #fasta = list(read_fasta(args.amplicon_file, args.minseqlen))
    jsp = list(dereplication_fulllength(args.amplicon_file, args.minseqlen, args.mincount))
    #jsp2 = abundance_greedy_clustering(args.amplicon_file, args.minseqlen, args.mincount, args.chunk_size, args.kmer_size)
    
    #print(list(jsp))
    
    #print(fasta)

#==============================================================
# Chimera removal section
#==============================================================

def get_unique(ids):
    return {}.fromkeys(ids).keys()

def common(lst1, lst2): 
    return list(set(lst1) & set(lst2))

def get_chunks(sequence, chunk_size):
    """Split sequences in a least 4 chunks
    """
    pass

def cut_kmer(sequence, kmer_size):
    """Cut sequence into kmers"""
    pass

def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    pass

def detect_chimera(perc_identity_matrix):
    pass

def search_mates(kmer_dict, sequence, kmer_size):
    pass

def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    pass


if __name__ == '__main__':
    main()
