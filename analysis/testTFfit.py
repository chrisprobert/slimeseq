#
#
#	file: testTFfit.py
#	date: 05/04/14
#	author: Chris Probert
#	
#	Purpose:
#	Test the fit of transcription factors on a set of upregulated genes
#	using a parameterized model of transcription factor binding likelihood.
#
#	Description:
#	Use a multinomial model of binding probabilities generated from an
# 	MLE of ChIP-seq tracks from many transcription factors to account
#	for position-specific bias. For each transcription factor being
#	tested, calculate a likelihood ratio based on the expected number
#	of binding sites within the threshold of the set of upregulated
#	genes being tested, versus the observed number of binding sites
#	within the threshold limit of DE genes for that TF.
#
#

import os, glob, sys
import numpy as np
import scipy.misc as scm

# path to the input data. All ".gff3" files are used.
PATH_TO_GFFS = "../input/modEncode/refinedPeaks/"
GFF_FILES = glob.glob(PATH_TO_GFFS + "*.gff3")

# path to wormbase genome file
PATH_TO_GENEFILE = "../genome/ws220bugfixGenes.txt"

# path to reference genome. Used for getting chromosome sizes.
PATH_TO_GENOME = "../genome/c_elegans.WS220.genomic.fa"

# path to directory containing differential expression tables
PATH_TO_DIFFEXPR_FILES = "../input/diffexpr/"
DE_FILES = glob.glob(PATH_TO_DIFFEXPR_FILES + "*.results.txt")

# window size for building null model of TF binding activity
WINDOW_SIZE = 200000

# binding range for how far a TF hit can be from a gene (in bp)
TF_RANGE = 2000

debug = True


def main() :

	print("getting chromosome lengths..")
	CHR_SIZES = getChrLengthsFromGenome(PATH_TO_GENOME)

	print("reading gene metadata..")
	genes = GeneCollection(PATH_TO_GENEFILE)
	if debug: print("Num genes total: %i" % len(genes.genes))

	print("reading gene DE data..")
	genes = updateGenesDEData(genes, DE_FILES)
	if debug: print("Num DE genes total: %i" % len([g for g in genes.genes if g.DE == True]))

	print("creating background probability model..")
	GFFs = [ GFFfile(filename) for filename in GFF_FILES ]
	allTFHits = []
	for GFF in GFFs :
		for hit in GFF.hits : allTFHits.append(hit)
	PM = probabilityModel(allTFHits, WINDOW_SIZE, CHR_SIZES, genes)

	print("calculating by-TF probabilities..")
	print("%s\t%s\t%s\t%s" % ("TF", "TotalHits", "Obs. DE-Hits", "Exp. DE-hits"))
	DEGenes = [g for g in genes.genes if g.DE == True]
	for GFF in GFFs :
		thisTFName = GFF.tf
		thisTFTotalHits = len(GFF.hits)
		thisTFOnTargetHits = []
		for hit in GFF.hits :
			onTarget = False
			for DEGene in DEGenes :
				if DEGene.chrm == hit.chr and np.abs(DEGene.start - hit.start) < TF_RANGE :
					onTarget = True
			if onTarget: thisTFOnTargetHits.append(hit)

		LikelihoodRatio = PM.calculateOverallLikelihoodRatio(thisTFOnTargetHits, thisTFTotalHits, DEGenes)
		print("%s\t%i\t%i\t%.8f" % (thisTFName, thisTFTotalHits, len(thisTFOnTargetHits), LikelihoodRatio[0]))

# class to hold contents of a GFF file for manipulation
class GFFfile :

	# ctor to create a GFF object from a given file path
	def __init__(self, filePath) :

		self.hits = []
		for line in open(filePath, "rU") :
			# check this line contains a valid tf hit entry
			if line[0] == "#" : continue
			if len(line.strip().split()) < 6 : continue
			self.hits.append(TFhit(line))
		self.tf = self.hits[0].tf

# class to hold a single TF hit from a GFF file
class TFhit :

	# ctor to take a line from a GFF file
	# example: III	EOR-1_L3_peaks	TF_binding_site	7524848	7525194	2.59E-08
	def __init__(self, line) :

		l = line.strip().split()

		self.chr = l[0]
		self.tf = l[1].split('_')[0]
		self.start = int(l[3])
		self.end = int(l[4])
		self.mid = int((self.end - self.start) / 2) + self.start
		self.pval = float(l[5])

# class to hold information on one gene
class Gene :

	# ctor to take information from gene file
	def __init__(self, chromosome, startPos, endPos, WBGeneID, publicName, desc="") :
		self.chrm = chromosome
		self.start = startPos
		self.end = endPos
		self.WBGeneID = WBGeneID
		self.publicName = publicName
		self.desc = desc
		self.logFC = 0
		self.PValue = 0
		self.FDR = 0
		self.DE = False

	# if this gene is differentially expressed, add the
	# information from differential expression analysis.
	def addDEInfo(self, logFC, PValue, FDR) :
		self.logFC = logFC
		self.PValue = PValue
		self.FDR = FDR
		self.DE = True

# class to hold a collection of upregulated genes
class GeneCollection :

	# read in the gene file, and store a list of the genes
	def __init__(self, pathToGeneFile) :
		"""
		expected schema for gene file (tsv) :
			Gene WB ID
			Gene Public Name
			Start (bp)
			End (bp)
			Chr Name
			Strand
			Regulation Summary (Regulated by)
			Regulation Result (Regulated by)
			Regulator Gene (WB Gene ID)
			Regulator Gene (Public Name)
			GO ID (merged)
		"""
		self.genes = []
		firstLine = False # keep track of whether we've parsed the header
		for line in open(pathToGeneFile, "rU") :

			if not firstLine and line[0] == "G" :
				firstLine = True
				continue
			try :
				l = line.split("\t")
				WBID = l[0]
				PubName = l[1]
				stPos = int(l[3])
				endPos = int(l[4])
				chrm = l[2]

				self.genes.append(Gene(chrm, stPos, endPos, WBID, PubName))

			except :
				continue
				#print("Error parsing gene file with line:\n%s\nStopping now." % line)

# class to hold genome-wide hit probability distribution
class probabilityModel :

	# ctor to initialize the probabilities based on a collection
	# of transcription factor hits
	def __init__(self, TFhits, window_size, genome_sizes, genes) :
		self.bins = {}
		self.probs = {}
		self.window_size = window_size
		self.genome_sizes = genome_sizes
		self.totalHits = 0 # keep track of how many total TF hits are included
		self.genes = genes

		for chrm in genome_sizes :
			if chrm == "MtDNA" : continue # skip mt DNA "chromosome"

			# calculate the number of bins based on the chromosome length
			chrmLength = genome_sizes[chrm]
			numBins = int(chrmLength / window_size) + 1

			self.bins[chrm] = [0 for i in range(numBins)]
			self.probs[chrm] = [0 for i in range(numBins)]

		for hit in TFhits :
			# use the midpoint of the hit to determine bin location :
			binNum = int(hit.mid / window_size)
			self.bins[hit.chr][binNum] += 1
			self.totalHits += 1

		for chrm in self.bins :
			for binNum in range(len(self.bins[chrm])) :
				self.probs[chrm][binNum] = float(self.bins[chrm][binNum]) / float(self.totalHits)

		#print self.probs["I"]
		totalP = float(0)
		for chrm in self.bins :
			for binNum in range(len(self.bins[chrm])) :
				totalP += self.probs[chrm][binNum]
		#print("%.12f" % totalP)

	def getProb(self, thisTFhit) :
		binNum = int(thisTFhit.mid / self.window_size)
		return self.probs[thisTFhit.chr][binNum]

	def calculateOverallLikelihoodRatio(self, OnTargetHits, TotalHitCount, DEGenes) :

		HitProb = float(0.0000001)
		for onTargetHit in OnTargetHits :
			# each hit is p(falling in TFRange | in window) * p(falling in window)
			HitProb += ((2*TF_RANGE) / WINDOW_SIZE) * self.getProb(onTargetHit)

		ExpHits = TotalHitCount * HitProb

		LikelihoodRatio = len(OnTargetHits) / ExpHits

		return [ExpHits, LikelihoodRatio]

# method to parse a genome fasta and get chromsome lengths
def getChrLengthsFromGenome(genomePath) :
	chr_sizes = {}
	curCount = 0
	curChr = ""
	for line in open(genomePath, "rU") :
		if line[0] == ">" :
			if curCount > 0 :
				chr_sizes[curChr] = curCount
			curCount = 0
			curChr = line.strip().split("_")[1]
		else :
			curCount += len(line.strip())
	chr_sizes[curChr] = curCount
	return chr_sizes

# update the GeneCollection with a list of DE text files
# return the updated GeneCollection
def updateGenesDEData(genes, DEfiles) :
	
	# expected format for DE file: WBGeneID	logConc	logFC	PValue	FDR

	for DEFileName in DEfiles :
		for line in open(DEFileName, "rU") :

			if line[0] == "l" : continue # skip header line

			l = line.strip().split()
			WBID = l[0].split(":")[0]
			logFC = float(l[2])
			PValue = float(l[3])
			FDR = float(l[4])

			isSet = False
			for i in range(len(genes.genes)) :
				if genes.genes[i].WBGeneID == WBID :
					if isSet : print(">> Found two gene objects with WBGeneID %s :-(" % WBID)
					isSet = True

					# only update this gene if it is not already marked as DE, or has a higher P-Value
					if genes.genes[i].DE == False or (genes.genes[i].PValue > PValue) :
						genes.genes[i].addDEInfo(logFC, PValue, FDR)

	return genes


if __name__ == "__main__" : main()