import sys
#sys.path.append("/home/athma/.local/lib/python3.6")
import pdb
import os
import subprocess
import argparse
import collections
from os.path import isfile, join
import gzip
import getopt
import re
import math
import random
import pysam
import pybedtools
import numpy as np
import scipy
from scipy import optimize
import pymc3 as pm
import pandas as pd

##### INITIAL PARSING + GETTING READS #####

#import dill
#dill.load_session('globalsave.pkl')

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

def getExons(bedname):
    if bedname[-2:] == "gz": bedfile = gzip.open(bedname)
    else: bedfile = open(bedname)
    beddict = collections.defaultdict(lambda:collections.defaultdict(dict))
    exoncount = 0
    while('TRUE'):
        bedline = bedfile.readline().strip().split('\t')
        if not bedline or bedline[0] == '': break
        gene = bedline[3].split(';')[1]
        exon = bedline[3].split(';')[0]
        exoncount += 1
        info = ';'.join(bedline[3].split(';')[2:7])
        start = bedline[1]
        end = bedline[2]
        strand = bedline[5]
        beddict[gene][exon]['info'] = info
        beddict[gene][exon]['start'] = start
        beddict[gene][exon]['end'] = end
        beddict[gene][exon]['strand'] = bedline[5]
    print("... " + str(exoncount) + " exons processed in " +str(len(beddict))+ " genes.")
    return(beddict)

def getJuncBed(bed, juncbam, readtype, readstrand):
    juncbam_name = juncbam[:-4]
    #print("basename: " +juncbam_name)
    #print("readtype: " +readtype)
    #print("readstrand: " +readstrand)
    if readstrand == 'fr-unstrand':
        # single bam file, no strand
        print("NOTE: unstranded reads")
        command = "intersectBed -abam " + juncbam + " -b " + bed + " -bed -wo -split | gzip > " + juncbam_name + "_exonjuncs.bed.gz"
        #command = "intersectBed -abam " + juncbam + " -b " + bed + " -bed -wo -split > " + juncbam_name + "_exonjuncs.bed"
        subprocess.call(command, shell=True)
    if readtype == 'single' and readstrand == 'fr-firststrand':
        # single bam file, read1=-S
        print("NOTE: fr-firststrand SE reads")
        command = "intersectBed -abam " + juncbam + " -b " + bed + " -bed -wo -S -split | gzip > " + juncbam_name + "_exonjuncs.bed.gz"
        subprocess.call(command, shell=True)
    if readtype == 'single' and readstrand == 'fr-secondstrand':
        # single bam file, read1=-s
        print("NOTE: fr-secondstrand SE reads")
        command = "intersectBed -abam " + juncbam + " -b " + bed + " -bed -wo -s -split | gzip > " + juncbam_name + "_exonjuncs.bed.gz"
        subprocess.call(command, shell=True)
    if readtype == 'paired' and readstrand == 'fr-firststrand':
        # two bams, read1=-S, read2=-s-
        print("NOTE: fr-firststrand PE reads")
        read1command = "intersectBed -abam " + juncbam_name + "_read1.bam -b " + bed + " -bed -wo -S -split > " + juncbam_name + "_read1.bed"
        read2command = "intersectBed -abam " + juncbam_name + "_read2.bam -b " + bed + " -bed -wo -s -split > " + juncbam_name + "_read2.bed"
        joincommand = "cat " + juncbam_name + "_read*.bed | gzip > " +juncbam_name +"_exonjuncs.bed.gz"  
        rmcommand = "rm -f " + juncbam_name + "_read*.bed"
        for command in (read1command, read2command, joincommand, rmcommand):
            subprocess.call(command, shell=True)
    if readtype == 'paired' and readstrand == 'fr-secondstrand':
        # two bams, read1=-s, read2=-S
        print("NOTE: fr-secondstrand PE reads")
        read1command = "intersectBed -abam " + juncbam_name + "_read1.bam -b " + bed + " -bed -wo -s -split > " + juncbam_name + "_read1.bed"
        read2command = "intersectBed -abam " + juncbam_name + "_read2.bam -b " + bed + " -bed -wo -S -split > " + juncbam_name + "_read2.bed"
        joincommand = "cat " + juncbam_name + "_read*.bed | gzip > " +juncbam_name +"_exonjuncs.bed.gz"  
        rmcommand = "rm -f " + juncbam_name + "_read*.bed"
        for command in (read1command, read2command, joincommand, rmcommand):
            subprocess.call(command, shell=True)

def getExonReads(juncbam, overlap):
    juncbam_name = juncbam[:-4]
    bedfile = gzip.open(juncbam_name + "_exonjuncs.bed.gz", 'rt')
    #bedfile = open(juncbam_name + "_exonjuncs.bed")
    readstartdict = collections.defaultdict(dict)
    readenddict = collections.defaultdict(dict)
    while('TRUE'):
        bedline = bedfile.readline().strip().split('\t')
        if not bedline or bedline[0] == '': break
        # check that read overlaps both ends of split region by at least overlap
        segments = bedline[10].split(',')
        if int(segments[0]) < int(overlap) or int(segments[1]) < int(overlap): continue
        # check that read overlaps this exon by at least overlap
        if int(bedline[18]) > (int(overlap) - 1):
            exon = bedline[15].split(';')[0] +';'+ bedline[15].split(';')[1]
            readstart = bedline[1]
            readend = bedline[2]
            if bedline[17] == '-':
                readstart = bedline[2]
                readend = bedline[1]
            # beddict[exon_name][read_name][read_position] = value
            readstartdict[exon][bedline[3]] = readstart
            readenddict[exon][bedline[3]] = readend
    return(readstartdict, readenddict)

##### HIT INDEX FUNCTION #####

def calculateMetric(beddict, startdict, enddict, readnum):
	#rows = []
	count = 0
	print("HIT for exons...",end='', flush=True)
	for gene in beddict:
		for exon in beddict[gene].keys():
			count += 1
			if count % 1000 == 0: print(str(count) +"..", end='', flush=True)
			# get data
			exonname = exon +';'+ gene
			#exstart = int(exon.split(':')[1].split('-')[0])
			#exend = int(exon.split(':')[1].split('-')[1])
			exstart = int(beddict[gene][exon]['start'])
			exend = int(beddict[gene][exon]['end'])
			# find reads starting/ending in exon
			covrange = range(exstart, exend)
			startreads = [k for k,v in startdict[exonname].items() if int(v) in covrange]
			endreads = [k for k,v, in enddict[exonname].items() if int(v) in covrange]
			# remove reads that start & end in exon
			bothstartend = set(startreads) & set(endreads)
			startreads_only = list(set(startreads) - bothstartend)
			endreads_only = list(set(endreads) - bothstartend)
			# calculate metric
			nright = len(startreads_only)
			nleft = len(endreads_only)
			HITindex = -2.0
			#if nright > 0 or nleft > 0:
			if nright + nleft > readnum:
				HITindex =  round(float(nleft - nright)/float(nright + nleft), 3)
			beddict[gene][exon]['name'] = exonname
			beddict[gene][exon]['nleft'] = nleft
			beddict[gene][exon]['nright'] = nright
			beddict[gene][exon]['HIT'] = HITindex
	print(str(count)+".")
	#df = pd.DataFrame(rows, columns=["name","gene","exon","nRIGHT","nTOTAL"])
	#return(beddict, df)
	return(beddict)

##### GENERATIVE MODEL FUNCTIONS #####

def classify_observations(N, D, k, omega):
	omega = omega[:-1]
	mode = np.array([1, 0.5, 0, 0])
	k = np.array(list(k) + [-1.0])
	k+=3
	omega = np.array(list(omega) + [1 - np.sum(omega)])

	alpha = mode * (k-2)+1
	beta = (1-mode)*(k-2)+1
	probs = []
	for i in range(len(k)):
		probs.append(omega[i] * scipy.stats.betabinom.pmf(np.array(D), np.array(N), alpha[i], beta[i]))
	
	probs = np.array(probs)
	probs = probs/probs.sum(0)
	return(probs)

def q_posterior(D, N, k_I, mesh=1000):
	q = np.linspace(0, 1, mesh)[1:-1]
	alpha = q * k_I
	beta = (1-q) * k_I
	lk = scipy.stats.betabinom.logpmf(D, N, alpha, beta)
	z = scipy.special.logsumexp(lk)
	return np.exp(lk-z)

def estimate_q(probs, D, N, k_I, mesh = 1000):
	summary_dict = {'post_mean': [],
					'prob_FI': [],
					'prob_IL': [],
					'CI_hyb_lo': [],
					'CI_hyb_hi': []}
	NH_means = np.array([1, 0.5, 0])
	center = int(mesh/2)-1
	for i in range(len(D)):
		post = q_posterior(D[i], N[i], k_I, mesh=1000)
		prob_FI = post[center:].sum() * probs[-1, i]
		prob_IL = post[:center].sum() * probs[-1, i]
		hybrid_mean = np.dot(post, np.linspace(0, 1, 1000)[1:-1])
		post_mean = (probs[:3,i] * NH_means).sum() + hybrid_mean*probs[-1,i]
		summary_dict['post_mean'].append(post_mean)
		summary_dict['prob_FI'].append(prob_FI)
		summary_dict['prob_IL'].append(prob_IL)
	for key in summary_dict.keys():
		summary_dict[key] = np.array(summary_dict[key])
	return summary_dict

def running_genmodel(exon_outname):
	df_exon = pd.read_csv(exon_outname, sep='\t')
	D = df_exon['nDOWN'].values
	N = (df_exon['nUP'] + df_exon['nDOWN']).values

	#with pm.Model() as model:
	#	# modes for each clas of metaexons held as constants
	#	mu = np.array([0, 0.5, 1])
	#	# prior reflecting how concentrated the beta distributions are
	#	# k at least 2
	#	k = pm.HalfNormal('k', 1000, shape=3)+3
	#	# convert mode and concentration parameters to alpha and beta parameters
	#	alpha = pm.Deterministic('alpha', mu*(k-2)+1)
	#	beta = pm.Deterministic('beta', (1-mu)*(k-2)+1)
	#	# omega is frequency of each class, modeled by uniform, Dirichlet prior
	#	omega = pm.Dirichlet('omega', np.array([1, 1, 1, 1]))
	#	# explictly define each distribution in the mixture
	#	betaF = pm.BetaBinomial.dist(alpha[2], beta[2], N)
	#	betaI = pm.BetaBinomial.dist(alpha[1], beta[1], N)
	#	betaL = pm.BetaBinomial.dist(alpha[0], beta[0], N)
	#	betaH = pm.BetaBinomial.dist(1, 1, N)
	#	# likelihood: observed number of downstream reads arises from mixture of beta-binom dists
	#	obs = pm.Mixture('D', w = omega, comp_dists = [betaF, betaI, betaL, betaH], observed = D)

	# for ADVI, minibatch fitting
	Dmini = pm.Minibatch(D, batch_size = 1000)
	Nmini = pm.Minibatch(N, batch_size = 1000)

	with pm.Model() as model_mini:
		# modes for each clas of metaexons held as constants
		mu = np.array([0, 0.5, 1])
		# prior reflecting how concentrated the beta distributions are
		# k at least 2
		k = pm.HalfNormal('k', 1000, shape=3)+3
		# convert mode and concentration parameters to alpha and beta parameters
		alpha = pm.Deterministic('alpha', mu*(k-2)+1)
		beta = pm.Deterministic('beta', (1-mu)*(k-2)+1)
		# omega is frequency of each class, modeled by uniform, Dirichlet prior
		omega = pm.Dirichlet('omega', np.array([1, 1, 1, 1]))
		# explictly define each distribution in the mixture
		betaF = pm.BetaBinomial.dist(alpha[2], beta[2], Nmini)
		betaI = pm.BetaBinomial.dist(alpha[1], beta[1], Nmini)
		betaL = pm.BetaBinomial.dist(alpha[0], beta[0], Nmini)
		betaH = pm.BetaBinomial.dist(1, 1, Nmini)
		# likelihood: observed number of downstream reads arises from mixture of beta-binom dists
		obs = pm.Mixture('D', w = omega, comp_dists = [betaF, betaI, betaL, betaH], observed = Dmini)

	with model_mini:
		mean_field = pm.fit(100000)

	mf_trace = mean_field.sample()
	# write back out lines with probabilities
	probs = []
	for i in range(mf_trace['k'].shape[0]):
		#if i%50 == 0: print(i)
		probs.append(classify_observations(N, D, mf_trace['k'][i,:], mf_trace['omega'][i,:]))
	probs = np.array(probs)
	probs = np.mean(probs, 0)
	# save to dataframe
	df_exon['PofF'] = np.round(probs[0], 5)
	df_exon['PofI'] = np.round(probs[1], 5)
	df_exon['PofL'] = np.round(probs[2], 5)
	df_exon['PofH'] = np.round(probs[3], 5)

	# getting probabilities of FI and IL
	mean_kI = mf_trace['k'][0].mean()
	summary = estimate_q(probs, D, N, mean_kI)
	# save to dataframe
	df_exon['PofFI'] = np.round(summary['prob_FI'], 5)
	df_exon['PofIL'] = np.round(summary['prob_IL'], 5)
	df_exon['downstream_fraction'] = np.round(summary['post_mean'], 5)
	df_exon['HIT_postmean'] = np.round(-1*(np.array(summary['post_mean'])*2-1), 5)

	df_exon.to_csv(exon_outname, sep='\t', index=False)

	return(df_exon)

##### DOWNSTREAM PROCESSING #####

def piecewise_linear(x, x0, y0, k1, k2): 
	return np.piecewise(x, [x < x0, x >= x0], [lambda x:k1*x + y0-k1*x0, lambda x:k2*x + y0-k2*x0])

def edge_flagging(HITdict):
	flaggingdata = []
	# for each exon, get distance from farthest exon with reads
	for gene in HITdict:
		# get all exon starts/ends with sum of reads >= readnum
		exonReads = [k for k in HITdict[gene].keys() if (HITdict[gene][k]['HIT'] > -2)]
		if len(exonReads) == 0: continue
		exonnames = [HITdict[gene][k]['name'] for k in exonReads] 
		exonleft = [HITdict[gene][k]['nleft'] for k in exonReads]
		exonright = [HITdict[gene][k]['nright'] for k in exonReads]
		exonratio = [exonright[i]/(exonright[i]+exonleft[i]) for i in range(len(exonReads))]
		strand = HITdict[gene][exonReads[0]]['strand']
		exstarts = [int(HITdict[gene][k]['start']) for k in exonReads]
		exends = [int(HITdict[gene][k]['end']) for k in exonReads]
		if strand == "+":
			startsite = min(exstarts)
			start_distances = [i - startsite for i in exstarts]
			endsite = max(exends)
			end_distances = [i - endsite for i in exends]
		if strand == "-":
			startsite = max(exends)
			start_distances = [startsite - i for i in exends]
			endsite = min(exstarts)
			end_distances = [endsite - i for i in exstarts]
		# add to dataframe
		data = {'name': exonnames, 'ratio': exonratio, 'startdistance': start_distances, 'enddistance': end_distances}
		df = pd.DataFrame(data=data)
		flaggingdata.append(df)
		# add distance to HITdict
		for i in range(len(exonnames)):
			HITdict[gene][exonReads[i]]['startdistance'] = start_distances[i]
			HITdict[gene][exonReads[i]]['enddistance'] = end_distances[i]

	# concatenate dataframes
	flaggingdata = pd.concat(flaggingdata, ignore_index=True)
	# piecewise regression to find deviation point (using genmodel N & D)
	ratio = np.array(flaggingdata['ratio'].values, dtype=float)
	distance = np.array(flaggingdata['startdistance'].values, dtype=float)
	optimresult = optimize.curve_fit(piecewise_linear, distance, ratio)
	param = optimresult[0][0]
	# flag exons smaller than deviation point (loop through HITdict)
	return(HITdict, param)

def writeMetrics(HITdict, param, outname):
    outexon = open(outname, 'w')
    outexon.write('exon\tgene\tstrand\tnTXPT\tnFE\tnINTERNAL\tnLE\tnSINGLE\tnUP\tnDOWN\tHITindex\tdist_to_TSS\tdist_to_TES\tedge\n')
    for gene in HITdict:
    	for exon in HITdict[gene].keys():
    		if HITdict[gene][exon]['HIT'] == -2: continue
    		# HITdict variables
    		infolist = [x.split(':')[1] for x in HITdict[gene][exon]['info'].split(';')]
    		exonname = HITdict[gene][exon]['name']
    		edge = 'no'
    		if HITdict[gene][exon]['startdistance'] <= param and HITdict[gene][exon]['HIT'] != -1 and HITdict[gene][exon]['HIT'] != 1:
    			edge = 'yes'
    		outexon.write(exon +'\t'+ gene +'\t'+ HITdict[gene][exon]['strand'] +'\t'+ '\t'.join(infolist) +'\t'+ 
    			str(HITdict[gene][exon]['nleft']) +'\t'+ str(HITdict[gene][exon]['nright']) +'\t'+ \
    			# HIT index
    			str(HITdict[gene][exon]['HIT']) +'\t'+ \
    			#str(HITdict[gene][exon]['bootpval'][0]) +'\t'+ str(HITdict[gene][exon]['bootpval'][1]) + \
    			# edge flag
    			str(HITdict[gene][exon]['startdistance']) +'\t'+ str(HITdict[gene][exon]['enddistance']) +'\t'+ edge +'\n')
    outexon.close()

##### IDENTIFY TERMINAL EXONS #####

def readMetrics(metrics):
	df = pd.read_csv(metrics, sep='\t')
	return(df)

def readParameters(parameters):
	paramfile = open(parameters)
	paramdict = dict()
	while('TRUE'):
		paramline = paramfile.readline().strip()
		if not paramline or paramline[0] == '': break
		if paramline[0] == "#": continue
		paramdict[paramline.split('\t')[0]] = paramline.split('\t')[1]
		# keys: HITterminal, HIThybrid, HITpval, HIT_CI, prob_med, prob_high
	return(paramdict)

def probabilistic_round(x):
    return int(math.floor(x + random.random()))

def sampleReads(allreads, nBOTH):
	nchoice = np.array(np.random.choice(allreads, nBOTH, replace=True))
	nchoice_up = np.sum(nchoice < 0)
	nchoice_down = np.sum(nchoice > 0)
	HIT = round(float(nchoice_up - nchoice_down) / float(nBOTH), 3)
	return(HIT)

def bootstrap(exon, n_boot, cutoff):
	nBOTH = int(exon.nUP) + int(exon.nDOWN)
	### bootstrap from null internal
	# create upstream "reads": negative numbers
	nUP_internal = probabilistic_round(nBOTH/2)
	up_internal = list(range((nUP_internal+1)*-1,0))
	# create downstream "reads": positive numbers
	nDOWN_internal = probabilistic_round(nBOTH/2)
	down_internal = list(range(1,(nDOWN_internal+2)))
	### bootstrap from confidence interval
	up_CI = list(range((int(exon.nUP)+1)*-1, 0))
	down_CI = list(range(1,(int(exon.nDOWN) + 2)))
	# all reads
	allreads_internal = up_internal + down_internal
	allreads_CI = up_CI + down_CI
	# get bootstrapped HITindices
	bootHITs_internal = []
	bootHITs_CI = []
	for boot in range(0, n_boot): 
		bootHITs_internal.append(sampleReads(allreads_internal, nBOTH))
		bootHITs_CI.append(sampleReads(allreads_CI, nBOTH))
	bootHITs_internal = np.array(bootHITs_internal)
	bootHITs_CI = np.array(bootHITs_CI)
	# calculate p-value - INTRON
	if float(exon.HITindex) == 0.0:
		pval_internal = np.sum(np.abs(bootHITs_internal) > 0.0) / float(n_boot)
	if float(exon.HITindex < 0.0):
		pval_internal = np.sum(bootHITs_internal <= exon.HITindex) / float(n_boot)
	if float(exon.HITindex) > 0.0:
		pval_internal = np.sum(bootHITs_internal >= exon.HITindex) / float(n_boot)
	# calculate CIs
	## calculate 75% confidence intervals
	CI75 = str(round(np.percentile(bootHITs_CI, 12.5), 3)) +','+ str(round(np.percentile(bootHITs_CI, 87.5), 3))
	## calculate 90% confidence intervals
	CI90 = str(round(np.percentile(bootHITs_CI, 5), 3)) +','+ str(round(np.percentile(bootHITs_CI, 95), 3))
	## calculate 95% confidence intervals
	CI95 = str(round(np.percentile(bootHITs_CI, 2.5), 3)) +','+ str(round(np.percentile(bootHITs_CI, 97.5), 3))
	# calculate p-value - CI
	pval_CI = np.sum(np.abs(bootHITs_CI) < abs(cutoff)) / float(n_boot)
	# return
	return CI75, CI90, CI95, pval_CI, pval_internal

def significance(HITcombo, n_boot, cutoff):
	HITcombo[['CI_75', 'CI_90', 'CI_95', 'pval_CI', 'pval_internal']] = HITcombo[['nUP', 'nDOWN', 'HITindex']].apply(bootstrap, axis=1, result_type="expand", n_boot=n_boot, cutoff=cutoff)
	#HITcombo['pval_IntronNull'] = HITcombo[['nUP', 'nDOWN', 'HITindex']].apply(bootIntronNull, axis=1, n_boot=n_boot)
	#HITcombo[['CI_75', 'CI_90', 'CI_95', 'pval_CI']] = HITcombo[['nUP','nDOWN']].apply(bootCI, axis=1, result_type="expand", n_boot=n_boot, cutoff=cutoff)
	HITcombo.to_csv(outname, sep='\t', index=False)	
	return(HITcombo)

def call_terminal(HITcombo, paramdict, outname):
	# set all to internal
	HITcombo['ID'] = "internal"
	# FI_med
	HITcombo.loc[(HITcombo.HITindex <= float(paramdict['HIThybrid'])*-1) & (HITcombo.PofFI >= float(paramdict['prob_med'])), 'ID'] = 'FirstInternal_medium'
	# FI_high
	HITcombo.loc[(HITcombo.ID == 'FirstInternal_medium') & (HITcombo.PofFI >= float(paramdict['prob_high'])), 'ID'] = 'FirstInternal_high'
	# IL_med
	HITcombo.loc[(HITcombo.HITindex >= float(paramdict['HIThybrid'])) & (HITcombo.PofIL >= float(paramdict['prob_med'])), 'ID'] = 'InternalLast_medium'
	# IL_high
	HITcombo.loc[(HITcombo.ID == 'InternalLast_medium') & (HITcombo.PofFI >= float(paramdict['prob_high'])), 'ID'] = 'InternalLast_high'
	# first
	HITcombo.loc[(HITcombo.HITindex <= float(paramdict['HITterminal'])*-1) & (HITcombo.pval_internal <= float(paramdict['HITpval'])), 'ID'] = 'first'
	# last
	HITcombo.loc[(HITcombo.HITindex >= float(paramdict['HITterminal'])) & (HITcombo.pval_internal <= float(paramdict['HITpval'])), 'ID'] = 'last'
	# confidence intervals
	if paramdict['HIT_CI'] != 'none':
		CIname = 'CI' + str(paramdict['HIT_CI'])
		HITcombo.loc[(HITcombo[CIname].split(',')[0] < 0) & (HITcombo[CIname].split(',')[1] > 0), 'ID'] = 'internal'
	# add new column with assignments based on position in gene
	HITcombo['ID_position'] = HITcombo['ID'] 
	# first if upstream most exon
	HITcombo.loc[(HITcombo.dist_to_TSS == 0), 'ID_position'] = 'first'
	# last if downstream most exon
	HITcombo.loc[(HITcombo.dist_to_TES == 0), 'ID_position'] = 'last'
	# write file to csv
	HITcombo.to_csv(outname, sep='\t', index=False)	
	return()

##### CALCULATE PSI #####

def calculatePSI(HITidentify, edge, outname):
	afepsidata = []
	alepsidata = []
	outAFE = open(outname +'.AFEPSI', 'w')
	#outAFE.write('exon\tgene\tstrand\tnFE\tnLEFT\tnRIGHT\tHITindex\tsum_R-L\tAFEPSI\n')
	outALE = open(outname +'.ALEPSI', 'w')
	#outALE.write('exon\tgene\tstrand\tnLE\tnLEFT\tnRIGHT\tHITindex\tsum_L-R\tALEPSI\n')
	# list of all genes
	genelist = HITidentify.gene.unique()
	for genehere in genelist:
		#print(genehere)
		dfhere = HITidentify.loc[HITidentify.gene == genehere]
		# get exons
		afehere = dfhere.loc[(dfhere.ID == 'first') | (dfhere.ID == 'FirstInternal_medium') | (dfhere.ID == 'FirstInternal_high')]
		if edge: afehere = afehere.loc[afehere.edge == 'no']
		alehere = dfhere.loc[(dfhere.ID == 'last') | (dfhere.ID == 'InternalLast_medium') | (dfhere.ID == 'InternalLast_high')]
		if edge: alehere = alehere.loc[alehere.edge == 'no']
		### CALCULATE AFE
		RLlist = afehere['nDOWN'] - afehere['nUP']
		if len(afehere) > 0 and sum(RLlist) > 0:
			if len(afehere) == 1: AFEPSI = [1.0]
			if len(afehere) > 1: AFEPSI = [float(x)/float(sum(RLlist)) for x in RLlist]
			# add to dataframe
			afehere['sumR-L'] = sum(RLlist)
			afehere['AFEPSI'] = AFEPSI
			# write afe psi values
			afehereparse = afehere[['gene', 'exon', 'strand', 'nTXPT','nFE', 'nUP', 'nDOWN', 'HITindex', 'sumR-L', 'AFEPSI']]
			afepsidata.append(afehereparse)
		### CALCULATE ALE
		LRlist = alehere['nUP'] - alehere['nDOWN']
		if len(alehere) > 0 and sum(LRlist) > 0:
			if len(alehere) == 1: ALEPSI = [1.0]
			if len(alehere) > 1: ALEPSI = [float(x)/float(sum(LRlist)) for x in LRlist]
			# add to dataframe
			alehere['sumL-R'] = sum(LRlist)
			alehere['ALEPSI'] = ALEPSI
			# write afe psi values
			alehereparse = alehere[['gene', 'exon', 'strand', 'nTXPT','nLE', 'nUP', 'nDOWN', 'HITindex', 'sumL-R', 'ALEPSI']]
			alepsidata.append(alehereparse)
	# write AFE
	afepsidata = pd.concat(afepsidata, ignore_index=True)
	afepsidata.to_csv(outAFE, sep='\t', index=False)	
	# write ALE
	alepsidata = pd.concat(alepsidata, ignore_index=True)
	alepsidata.to_csv(outALE, sep='\t', index=False)	

if __name__ == '__main__':
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='Classify and quantify exons from a bam file.', add_help = True)
	parser.add_argument('--junctionReads', action='store_true', default=False, help='Extract junction reads', required = False)
	parser.add_argument('--HITindex', action='store_true', default=False, help='Calculate HITindex', required = False)
	parser.add_argument('--classify', action = 'store_true', help='Classify terminal, hybrid, and internal exons', required = False)
	parser.add_argument('--calculatePSI', action = 'store_true', help='Calculate PSI values', required = False)
	parser.add_argument('--outname', type=str, metavar='', help='name of file(s) for final metric. required for everything except --junctionReads.', required=False)
	### read information
	group_reads = parser.add_argument_group('read information')
	group_reads.add_argument('--bam', type = str, metavar='', help = 'original bam from which to extract junction reads. required if --junctionReads', required=False, default="None")
	group_reads.add_argument('--juncbam', type = str, metavar='', help = 'junction read bam. required if --junctionReads or --HITindex', required = False, default = "None")
	group_reads.add_argument('--readtype', type=str, help = 'type of read', default='paired', required=False, choices=['single','paired'])
	group_reads.add_argument('--readstrand', type=str, help = 'directionality of RNA-seq data', default='fr-firststrand', required=False, choices=['fr-unstrand','fr-firststrand','fr-secondstrand'])
	### exon information
	group_exons = parser.add_argument_group('exon information')
	group_exons.add_argument('--bed', type = str, metavar='', help = 'bed file with merged/annotated exons. Output from HITindex_annotate.py. required if --HITindex', required = False, default = "None")
	## HITindex parameters
	group_param = parser.add_argument_group('HITindex', 'parameters for running HIT index')
	group_param.add_argument('--overlap', type = int, metavar='', help = 'overlap of split read with exon region (nt)', required = False, default = 10)
	group_param.add_argument('--readnum', type = int, metavar='', help = 'minimum number of reads for confidence in HITindex (sum of R + L)', required = False, default = 2)
	### thresholds for calling terminal exons
	group_term = parser.add_argument_group('classify', 'information for classifying exon types') 
	group_term.add_argument('--metrics', metavar='', type = str, help = 'HITindex output file, required if --HITindex is not specified.', required = False, default = "None")
	group_term.add_argument('--parameters', metavar='', type = str, help = 'file specifying HITindex and generative model thresholds for classifying exons.', required=False, default="HIT_identity_parameters.txt")
	group_term.add_argument('--bootstrap', type = int, metavar='', help = 'bootstrapping iterations to get confidence intervals and p-values', required = False, default = 1000)
	### PSI values
	group_psi = parser.add_argument_group('psi', 'parameters for calling PSI values')
	group_psi.add_argument('--metricsID', type = str, metavar='', help = 'HITindex identification output file, required if --classify is not specified.', required = False, default = "None")
	group_psi.add_argument('--edge', type = bool, metavar='', help = 'exclude exons flagged as being affected by the edge effect from PSI calculations', required = False, default = False)
	# confidence information
	args = parser.parse_args()

	if not args.junctionReads and not args.HITindex and not args.classify and not args.calculatePSI:
		sys.exit("ERROR! Need to include at least one function: --junctionReads, --HITindex, --classify, and/or --calculatePSI")

    #pdb.set_trace()
	# dill.dump_session('globalsave.pkl')

	if args.junctionReads:
		if args.bam == "None" or args.juncbam == "None":
			sys.exit("ERROR! Need to include --bam and --juncbam when running --junctionReads")
		print('Obtaining junctions...')
		if args.readtype == 'single' or args.readstrand == 'fr-unstrand':
			JunctionBam_SE(args.bam, args.juncbam)
		if args.readtype == 'paired' and args.readstrand != 'fr-unstrand':
			JunctionBam_PE(args.bam, args.juncbam)
		print('... junctions obtained.')

	if args.HITindex:
		if args.juncbam == "None" or args.bed == "None":
			sys.exit("ERROR! Need to include --juncbam and --bed when running --HITindex")
		print("Calculating metric...")
		beddict = getExons(args.bed)
		print("... got exons...")
		getJuncBed(args.bed, args.juncbam, args.readtype, args.readstrand)
		print("... got junctions...")
		startdict, enddict = getExonReads(args.juncbam, args.overlap)
		print("... got counts...")
		HITdict = calculateMetric(beddict, startdict, enddict, args.readnum)
		print("... HITindex calculated.")
		HITdict, param = edge_flagging(HITdict)
		print("... edge effect exons flagged.")
		exon_outname = args.outname +'.exon'
		writeMetrics(HITdict, param, exon_outname)
		print("... metrics written.")
		genmodel = running_genmodel(exon_outname)
		#genmodel_dict = genmodel.set_index('name').T.to_dict('list')
		# 'name': ['gene, 'exon', 'nRIGHT', 'nTOTAL', 'PofF', 'PofI', 'PofL', 'PofH', 'PofFI', 'PofIL', 'downstream_fraction', 'HIT_postmean']
		print("... generative model calculated.")

	if args.classify:
		if not args.HITindex:
			if args.metrics == "None":
				sys.exit("ERROR! Need to include --metrics specifying HITindex output file to use")
			exon_outname = args.metrics
		print("Classifying exons...")
		HITcombo = readMetrics(exon_outname)
		print("...read exons.")
		paramdict = readParameters(args.parameters)
		print("...read parameters.")
		HITcombo = significance(HITcombo, args.bootstrap, float(paramdict['HIThybrid']))
		print("...calculated significance.")
		HITidentify = call_terminal(HITcombo, paramdict, exon_outname)
		print("... exons identified.")

	if args.calculatePSI:
		if not args.classify:
			if args.metricsID == "None":
				sys.exit("ERROR! Need to include --metricsID specifying HITindex output file to use with exons identified")
			HITidentify = readMetrics(args.metricsID)
		print("Calculating PSI...")
		calculatePSI(HITidentify, args.edge, args.outname)
		print("... PSI calculated.")
