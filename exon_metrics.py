import sys
import argparse
import pdb
import os
import subprocess
import gzip

from exon_annot import getGenes, annotateGenes
from calculate_metric import getExons, getExonReads, calculateMetric, calculatePSI

def JunctionBam_SE(bam, outfile):
    bamname = bam[:-4]
    cmd_header = "samtools view -H " + bam + " > " + outfile + "_header.txt"
    cmd_junc = "samtools view -F 256 " + bam + " | awk '{if ($6 ~/N/) {print $0}}' | cat " + outfile + "_header.txt - | samtools view -bS - | samtools sort - -T " \
               + outfile + " -o " + outfile
    cmd_index = "samtools index " +outfile
    for command in (cmd_header, cmd_junc, cmd_index):
        subprocess.call(command, shell=True)

def JunctionBam_PE(bam, outfile):
    bamname = bam[:-4]
    cmd_header = "samtools view -H " + bam + " > " + outfile + "_header.txt"
    # read 1
    outfile1 = outfile[:-4] +"_read1.bam"
    cmd_junc = "samtools view -f 64 -F 256 " + bam + " | awk '{if ($6 ~/N/) {print $0}}' | cat " + outfile + "_header.txt - | samtools view -bS - | samtools sort - -T " \
               + outfile1 + " -o " + outfile1
    cmd_index = "samtools index " +outfile1
    for command in (cmd_header, cmd_junc, cmd_index):
        subprocess.call(command, shell=True)
    # read 2
    outfile2 = outfile[:-4] +"_read2.bam"
    cmd_junc = "samtools view -f 128 -F 256 " + bam + " | awk '{if ($6 ~/N/) {print $0}}' | cat " + outfile + "_header.txt - | samtools view -bS - | samtools sort - -T " \
               + outfile2 + " -o " + outfile2
    cmd_index = "samtools index " +outfile2
    for command in (cmd_junc, cmd_index):
        subprocess.call(command, shell=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # get junctions
    parser.add_argument('--junctions', help = 'Get junction reads for any bam file', action = 'store_true', required = False)
    parser.add_argument('--bam', type = str, help = 'bam from which to extract junction reads. required if --junctions', required = False, default = "NotSpec")
    parser.add_argument('--outfile_junc', type = str, help = 'name for output bam with junction reads. required if --junctions', required = False, default = "NotSpec")
    parser.add_argument('--readtype', type=str, help = 'type of read (default = paired)', default='paired', required=False, choices=['single','paired'])
    parser.add_argument('--readstrand', type=str, help = 'directionality of RNA-seq data (default = fr-firststrand)', default='fr-firststrand', required=False, choices=['fr-unstrand','fr-firststrand','fr-secondstrand'])
    # get metrics
    parser.add_argument('--metric', help = 'Calculate metric for all annotated exons', action = 'store_true', required = False)
    parser.add_argument('--bed', type = str, help = 'bed file with merged/annotated exons. same as --outfile_bed. required if --metric', required = False, default = "NotSpec")
    parser.add_argument('--juncbam', type = str, help = 'bam file with junction reads. same as --outfile_bam. required if --metric', required = False, default = "NotSpec")
    parser.add_argument('--overlap', type = int, help = 'overlap of split read with exon region (nt). default = 10.', required = False, default = 10)
    parser.add_argument('--readnum', type=int, help='minimum number of reads for confidence in HITindex (sum of R + L). default = 2', required = False, default = 2)
    parser.add_argument('--bootstrap', type = str, help = 'bootstrapping iterations to get p-value for metric confidence (within 0.1). default = 100', required = False, default = 100)
    parser.add_argument('--AFEcutoff', type=float, help='HITindex cutoff to call an exon a first exon. default = 0.2', required = False, default = 0.2)
    parser.add_argument('--ALEcutoff', type=float, help='HITindex cutoff to call an exon a last exon. defalut = -0.2', required = False, default = -0.2)
    parser.add_argument('--outname', type = str, help = 'name of file(s) for final metrics. required if --metric', required = False, default = "NotSpec")
    args = parser.parse_args()

    # get junctions
    if args.junctions == True:
        if args.bam != "NotSpec" and args.outfile_junc != "NotSpec":
            print('Obtaining junctions...')
            if args.readtype == 'single' or args.readstrand == 'fr-unstrand':
                JunctionBam_SE(args.bam, args.outfile_junc)
            if args.readtype == 'paired' and args.readstrand != 'fr-unstrand':
                JunctionBam_PE(args.bam, args.outfile_junc)
            print('... junctions obtained.')
        else:
            print("ERROR: when using --junctions, --bam AND --outfile_bam must be specified")

    # get metrics + PSI
    if args.metric == True:
        if args.bed !="NotSpec" and (args.juncbam !="NotSpec" or args.outfile_junc !="NotSpec") and args.outname !="NotSpec":
            print('Calculating metric...')
            if args.juncbam == "NotSpec": args.juncbam = args.outfile_junc
            beddict = getExons(args.bed)
            getJuncBed(args.bed, args.juncbam, args.readtype, args.readstrand)
            startdict, enddict = getExonReads(args.juncbam, args.overlap)
            HITdict = calculateMetric(beddict, startdict, enddict, args.AFEcutoff, args.ALEcutoff, args.outname)
            print('... metric calculated.')
            print('Calculating PSI...')
            calculatePSI(HITdict, args.AFEcutoff, args.ALEcutoff, args.readnum, args.outname)
            print('... PSI calculated.')
        else:
            print("ERROR: when using --metric, --bed AND --juncbam AND --outname must be specified")