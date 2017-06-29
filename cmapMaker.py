#!/usr/bin/env python2

# Modules need: 'biopython', 'pathos'
# Written by: Yuxuan (Andy) Yuan
# Email: yuxuan.yuan@research.uwa.edu.au
# Created date: 20/02/2017
# Last modification date: 28/06/2017
## things need to do: 1) creat a function to add customized enzymes

#Usage: python cmapMaker.py -f test.fa -e BssSI BspQI -l 0 -s 0 -t 2 -o . 

#======================= Modules ======================

from Bio import SeqIO
from Bio.Seq import Seq
import pathos.multiprocessing as mp
import argparse
import os
import sys
import re

#======================== Function ==========================

class cmapMaker (object):

    def __init__(self, fasta_file, enzs, min_len, min_nsite, threads, path):
        """
        This script is used to digest NGS sequences into cmap files required by BioNano analysis.
        """
        self.fasta_file = fasta_file
        self.name = fasta_file.rsplit('.',1)[0].split('/')[-1]
        self.enzs = enzs
        self.min_len = min_len
        self.min_nsite = min_nsite
        self.threads = threads
        self.path = path

    def check_fa(self):
        self.fasta_file = os.path.abspath(self.fasta_file)
        assert os.path.isfile(self.fasta_file) and os.access(self.fasta_file,os.R_OK),\
               "The fasta file {} doesn't exist or isn't readable".format(self.fasta_file)
    
    def check_path(self):
        self.path = os.path.abspath(self.path)
        if not os.path.isdir(self.path):
            print >> sys.stderr, '\nOops! It seems the path to output directory is not existent. Please check!\n'
            sys.exit()
    def checkPathW(self):
        if os.access(self.path, os.W_OK) is not True:
            print >> sys.stderr, '\nOops! It seems the outDir is not writable. Please check!.\n'
    

    def cmap_maker(self,enz):
        index = 0
        enzymes = {'BspQI':'GCTCTTC',
                    'BbvCI':'CCTCAGC',
                    'BsmI':'GAATGC',
                    'BsrDI':'GCAATG',
                    'bseCI':'ATCGAT',
                    'BssSI':'CACGAG',
                    'HindIII':'AAGCTT'} # users can modify this to add new enzymes using same format
        try:
            cmap_file='%s/%s_%s_%sKb_%slabels.cmap'%(self.path,self.name,enz,self.min_len,self.min_nsite)
            forwards = enzymes[enz]
            reverse = str(Seq(forwards).reverse_complement())
            with open (cmap_file,'w') as ref_cmap:
                ref_cmap.write('# CMAP File Version:\t0.1\n')
                ref_cmap.write('# Label Channels:\t1\n')
                ref_cmap.write('# Nickase Recognition Site 1:\t%s\n'%forwards)
                if enz !='BspQI' and enz !='HindIII':
                    ref_cmap.write('# Enzyme1:\tNb.%s\n'%enz)
                else:
                    if enz !='HindIII':
                        ref_cmap.write('# Enzyme1:\tNt.%s\n'%enz)
                    else:
                        ref_cmap.write('# Enzyme1:\t%s\n'%enz)
                ref_cmap.write('# Number of Consensus Nanomaps:\tTBD\n')
                ref_cmap.write('#h CMapId\tContigLength\tNumSites\tSiteID\tLabelChannel\tPosition\tStdDev\tCoverage\tOccurrence\n')
                ref_cmap.write('#f int\tfloat\tint\tint\tint\tfloat\tfloat\tint\tint\n')
                for seqs in SeqIO.parse(self.fasta_file,'fasta'):
                    seq = str(seqs.seq.upper())
                    seq_len = len(seq)
                    index+=1
                    if seq_len >= self.min_len*1000:
                        nsites = len(re.findall('%s|%s'%(forwards,reverse),seq))
                        if nsites >=self.min_nsite:
                            j=1
                            for o in re.finditer('%s|%s'%(forwards,reverse),seq):
                                ref_cmap.write('%s\t%.1f\t%d\t%d\t1\t%.1f\t1.0\t1\t1\n'%(index,seq_len,nsites,j,o.start()+1))
                                j+=1
                            ref_cmap.write('%s\t%.1f\t%d\t%d\t0\t%.1f\t0.0\t1\t0\n'%(index,seq_len,nsites,j,seq_len))
        except:
            print >> sys.stderr, '\nPlease check your input fasta file or the writability of your output directory!\n'
            os.remove('%s/%s_%s_%sKb_%slabels.cmap'%(self.path,self.name,enz,self.min_len,self.min_nsite))
            sys.exit()
        return
    def digestion(self):
        pool= mp.ProcessingPool(self.threads)
        enzs = self.enzs
        pool.map(self.cmap_maker, enzs)

if __name__=="__main__":
    """
    Parse the function and usage.
    """
    parser = argparse.ArgumentParser(description='Digest a fasta format file into cmap files using selected enzymes.')
    parser.add_argument('-v', '--version', action='version', version='1.01')
    parser.add_argument('-f', dest='fasta', help='a fasta format file that contains all genome sequences', type = str)
    parser.add_argument('-e', dest='enzyme', nargs='*', help='enzyme(s) to digest the NGS sequences [BspQI|BbvCI|BsmI|BsrDI|bseCI|BssSI|HindIII]. Can be one or multiple.', type = str) # users can modify this to add new enzymes using same format
    parser.add_argument('-l', dest='length', help='minimum length of the sequences (kb). Default: 20.', default = 20, type=int)
    parser.add_argument('-s', dest='site', help='minimum label of enzyme nicking sits on the sequences. Default: 5.', default = 5, type=int)
    parser.add_argument('-t', dest='threads', help='number of threads. Default: 1', default = 1, type=int)
    parser.add_argument('-o', dest='output', help='the output directory', type=str)
    args = parser.parse_args()
    if None not in [args.fasta, args.enzyme, args.output]:
        for enz in args.enzyme:
            if enz not in ['BspQI', 'BbvCI', 'BsmI', 'BsrDI', 'bseCI', 'BssSI', 'HindIII']: # users can modify this to add new enzymes using same format
                print >> sys.stderr, '\nPlease input enzyme(s) from [BspQI|BbvCI|BsmI|BsrDI|bseCI|BssSI|HindIII] with a whitespace between each other.\n' # users can modify this to add new enzymes using same format
                sys.exit()
            
        run=cmapMaker(args.fasta, args.enzyme, args.length, args.site, args.threads, args.output)
        run.check_fa()
        run.check_path()
        run.checkPathW()
        run.digestion()
                    
    else:
        print
        parser.print_help()
        print
        sys.exit()
