import sys
import pdb
import os
import subprocess
import argparse
import collections
from os.path import isfile, join
import gzip
import getopt
import re
import pysam
import pybedtools

# read in genes and annotations
def getGenes(gtf):
    if gtf[-2:] == "gz": gtffile = gzip.open(gtf)
    else: gtffile = open(gtf)

    # dictionary to store genes
    # structure: genedict[gene_name][transcript_name][exon_#][exon_coords]
    genedict = collections.defaultdict(lambda:collections.defaultdict(dict))

    while('TRUE'):
        gtfrow = gtffile.readline()
        if not gtfrow: break
        gtfline = gtfrow.split()
        if gtfline[0][:1] == "#": continue
        if gtfline[2] == "gene":
            gene = gtfline[9].split('\"')[1] +":"+ gtfline[0] +":"+ gtfline[3] +"-"+ gtfline[4] +":"+ gtfline[6]
        if gtfline[2] == "transcript":
            transcript = gtfline[13].split('\"')[1] +":"+ gtfline[0] +":"+ gtfline[3] +":"+ gtfline[4]
            ind = 0
        if gtfline[2] == "exon":
            genedict[gene][transcript][ind] = gtfline[3] +"_"+ gtfline[4]
            ind = ind+1
    print(str(len(genedict)) +" genes input.")
    return(genedict)

# merge exons for each gene, annotate exons with as terminal/internal, write out bed file
def annotateGenes(genedict, reverse, outfile, ss3buffer, ss5buffer):
    outfh = open(outfile, 'w')
    outfhbuffer = open(outfile +'_ss3-'+ str(ss3buffer) +'ss5-'+ str(ss5buffer) +'.buffer', 'w') 
    outfhconst = open(outfile +'_constituents', 'w')
    # list of possible exon types - terminal, internal, and singleton
    exontypes = ['FE', 'internal', 'LE', 'singleexon']
    
    for gene in genedict:
        print(gene)
        # gene info
        genesymbol = gene.split(":")[0]
        exons = []
        exons_const = []
        chr = gene.split(":")[1]
        strand = gene.split(":")[3]
        ntxpt = 'TXPT:' + str(len(genedict[gene]))
    
        for transcript in genedict[gene]:
            #get all exons within transcript
            exonshere = genedict[gene][transcript]
            nexonshere = len(exonshere)
            if nexonshere == 1: exnames = ['singleexon']
            if nexonshere > 1:
                if strand == '+': exnames = ['FE'] + ['internal']*(nexonshere - 2) + ['LE']
                if strand == '-': exnames = ['LE'] + ['internal']*(nexonshere - 2) + ['FE']
            if nexonshere > 1 and reverse == True:
                exnames = ['FE'] + ['internal']*(nexonshere - 2) + ['LE']
            for nex in range(0, nexonshere):
                ex = exonshere[nex].split('_')
                exons.append(chr +' '+ ex[0] +' '+ ex[1] +' '+ exnames[nex] +' . '+ strand)
                exons_const.append(chr +' '+ ex[0] +' '+ ex[1] +' '+ ex[0] +'-'+ ex[1]+ ' . '+ strand)
        
        # merge overlapping exons into single exon - maintaining information about exon type in original transcripts 
        bedscratch = pybedtools.BedTool('\n'.join(exons), from_string=True)
        bedsort = bedscratch.sort()
        bedsortmerge = bedsort.merge(s=True, c=4, o='collapse')
        bedlist = str(bedsortmerge).split('\n')[:-1]

        # same for constituent exon list
        bedscratch_const = pybedtools.BedTool('\n'.join(exons_const), from_string=True)
        bedsort_const = bedscratch_const.sort()
        bedsortmerge_const = bedsort_const.merge(s=True, c=4, o='collapse')
        bedlist_const = str(bedsortmerge_const).split('\n')[:-1]

        # write out bed of merged exons
        for ex in bedlist:
            exhere = ex.split('\t')
            exherename = exhere[0] +':'+ exhere[1] +'-'+ exhere[2]
            extypes = exhere[3].split(',')
            ntypeset = []
            # count types of exons
            for x in exontypes:
                ntypeset.append(x +':'+ str(extypes.count(x)))
           # get buffered regions by strand
            if strand == '+':
                exstart = int(exhere[1]) - ss3buffer
                exend = int(exhere[2]) + ss5buffer
            if strand == '-':
                exstart = int(exhere[1]) - ss5buffer
                exend = int(exhere[2]) + ss3buffer
            outfh.write(chr +'\t'+ exhere[1] +'\t'+ exhere[2] +'\t'+ exherename +';'+ genesymbol +';'+ ntxpt +';'+ ';'.join(ntypeset) +'\t0\t'+ strand +'\n')
            outfhbuffer.write(chr +'\t'+ str(exstart) +'\t'+ str(exend) +'\t'+ exherename +';'+ genesymbol +';'+ ntxpt +';'+ ';'.join(ntypeset) +'\t0\t'+ strand +'\n')

        # write out bed of merged exons
        for ex in bedlist_const:
            exhere = ex.split('\t')
            exherename = exhere[0] +':'+ exhere[1] +'-'+ exhere[2]
            exactual = exhere[3]
            outfhconst.write(chr +'\t'+ exhere[1] +'\t'+ exhere[2] +'\t'+ exherename +';'+ genesymbol +';'+ ntxpt +';'+ exactual +'\t0\t'+ strand +'\n')
            
    outfh.close()
    outfhbuffer.close()
    outfhconst.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='Identify and annotate metaexons from a gtf annotation set.')
    # inputs
    group_input = parser.add_argument_group('Input')
    group_input.add_argument('--gtf', dest = 'gtf', type = str, metavar = 'gtf', help = 'gtf to be indexed', required=True)
    # parameters
    group_param = parser.add_argument_group('Parameters')
    group_param.add_argument('--reverse', action="store_true", help = 'use if exons are sorted by transcriptional direction rather than reference coordinate', required = False)
    group_param.add_argument('--ss3buffer', type = int, metavar = '', help = 'intronic buffer region included upstream of 3ss of exon for counting reads. suggested = 50nt.', required = False, default = 0)
    group_param.add_argument('--ss5buffer', type = int, metavar = '', help = 'intronic buffer region included downstream of 5ss of exon for counting reads. suggested = 20nt.', required = False, default = 0)
    # outputs
    group_output = parser.add_argument_group('Output')
    group_output.add_argument('--outfile', type = str, metavar = 'output', help = 'name for output bed with merged/annotated exons', required = True)

    args = parser.parse_args()

    print('Indexing gtf...')
    genedict = getGenes(args.gtf)
    annotateGenes(genedict, args.reverse, args.outfile, args.ss3buffer, args.ss5buffer)
    print('... gtf indexed.')
