#!/usr/bin/env python2

# Written by: Andy Yuan
# Email: yuxuan.yuan@outlook.com
# Last modification date: 19/04/2017
## things need to do: 1) creat a function to add customized enzymes

from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import os
import sys
import re


def cmap_maker(fasta_file, enz, min_len, min_nsite, path):
    """
    This script is used to digest NGS sequences into a cmap file required by BioNano.
    Currently, six enzymes are provided, which are BspQI, BbvCI, Bsml, BsrDI, bseCI and BssSI.
    """
    name = fasta_file.rsplit('.',1)[0].split('/')[-1]
    fasta_file = os.path.abspath(fasta_file)
    path = os.path.abspath(path)
    if not os.path.isdir(path):
        print >> sys.stderr, '\nSomething is wrong with your output directory! Please check!'
        sys.exit()
    index = 0
    enzymes = {'BspQI':'GCTCTTC',
                'BbvCI':'CCTCAGC',
                'BsmI':'GAATGC',
                'BsrDI':'GCAATG',
                'bseCI':'ATCGAT',
                'BssSI':'CACGAG',
                'HindIII':'AAGCTT'} # users can modify this to add new enzymes using same format
    try:
        cmap_file='%s/%s_%s_%sKb_%slabels.cmap'%(path,name,enz,min_len,min_nsite)
        forwards = enzymes[enz]
        reverse = str(Seq(forwards).reverse_complement())
        with open (cmap_file,'w') as ref_cmap:
            ref_cmap.write('# CMAP File Version:\t0.1\n')
            ref_cmap.write('# Label Channels:\t1\n')
            ref_cmap.write('# Nickase Recognition Site 1:\t%s\n'%forwards)
            if enz !='BspQI':
                ref_cmap.write('# Enzyme1:\tNb.%s\n'%enz)
            else:
                ref_cmap.write('# Enzyme1:\tNt.%s\n'%enz)
            ref_cmap.write('# Number of Consensus Nanomaps:\tTBD\n')
            ref_cmap.write('#h CMapId\tContigLength\tNumSites\tSiteID\tLabelChannel\tPosition\tStdDev\tCoverage\tOccurrence\n')
            ref_cmap.write('#f int\tfloat\tint\tint\tint\tfloat\tfloat\tint\tint\n')
            for seqs in SeqIO.parse(fasta_file,'fasta'):
       	        seq = str(seqs.seq.upper())
       	        seq_len = len(seq)
       	        index+=1
       	        if seq_len >= min_len*1000:
                    nsites = len(re.findall('%s|%s'%(forwards,reverse),seq))
                    if nsites >=min_nsite:
                        j=1
                        for o in re.finditer('%s|%s'%(forwards,reverse),seq):
                            ref_cmap.write('%s\t%.1f\t%d\t%d\t1\t%.1f\t1.0\t1\t1\n'%(index,seq_len,nsites,j,o.start()+1))
                            j+=1
                        ref_cmap.write('%s\t%.1f\t%d\t%d\t0\t%.1f\t0.0\t1\t0\n'%(index,seq_len,nsites,j,seq_len))
    except:
        print >> sys.stderr, '\nPlease check your input fasta file or the writability of your output directory!\n'
        os.remove('%s/%s_%s_%sKb_%slabels.cmap'%(path,name,enz,min_len,min_nsite))
        sys.exit()



if __name__=="__main__":
    """
    Parse the function and usage.
    """
    parser = argparse.ArgumentParser(description='Digest a fasta format file into a cmap file using a selected enzyme.')
    parser.add_argument('-v', '--version', action='version', version='1.0')
    parser.add_argument('-f', dest='fasta', help='a fasta format file that contains all genome sequences', type = str)
    parser.add_argument('-e', dest='enzyme', help='an enzyme to digest the NGS sequences [BspQI|BbvCI|BsmI|BsrDI|bseCI|BssSI|HindIII]', type = str) # users can modify this to add new enzymes using same format
    parser.add_argument('-l', dest='length', help='minimum length of the sequences (kb). Default: 20.', default = 20, type=int)
    parser.add_argument('-s', dest='site', help='minimum label of enzyme nicking sits on the sequences. Default: 5.', default = 5, type=int)
    parser.add_argument('-o', dest='output', help='the output directory', type=str)
    args = parser.parse_args()
    if None not in [args.fasta, args.enzyme, args.output]:
        if args.enzyme in ['BspQI', 'BbvCI', 'Bsml', 'BsrDI', 'bseCI', 'BssSI', 'HindIII']: # users can modify this to add new enzymes using same format
            try:
                cmap_maker(args.fasta, args.enzyme, args.length, args.site, args.output)
            except:
                print >> sys.stderr, '\nSomething is wrong with you input, please check!\n'
        else:
            print >> sys.stderr, '\nPlease input an enzyme from [BspQI|BbvCI|BsmI|BsrDI|bseCI|BssSI|HindIII].\n' # users can modify this to add new enzymes using same format
            sys.exit()
    else:
        print
        parser.print_help()
        print
        sys.exit()
