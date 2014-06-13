import numpy as np

# class to hold contents of a GFF file for manipulation
class GFFfile :

	# ctor to create a GFF object from a given file path
	def __init__(self, filePath) :

		self.filePath = filePath
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
	def __init__(self, chromosome, strand, startPos, endPos, WBGeneID, publicName, desc="") :
		self.chrm = chromosome
		self.strand = strand
		self.start = startPos
		self.end = endPos
		self.WBGeneID = WBGeneID
		self.publicName = publicName
		self.desc = desc
		self.logFC = 0
		self.PValue = 1
		self.FDR = 1
		self.DE = False
		self.GroveOverlap = 0

	# if this gene is differentially expressed, add the
	# information from differential expression analysis.
	def addDEInfo(self, logFC, PValue, FDR) :
		self.logFC = logFC
		self.PValue = PValue
		self.FDR = FDR
		self.DE = True

	def addGroveOverlap(self) :
		self.GroveOverlap = 1

	# check if a given TF hit matches this gene
	def checkTFHitValid(self, TFHit, TF_RANGE_BEFORE, TF_RANGE_AFTER) :

		# sanity checks on the TF hit
		if TFHit.mid <= 0 : return False
		if not TFHit.mid : return False

		if TFHit.chr != self.chrm : return False

		# the bounds to compare TF hit position against
		leftBound = 0
		rightBound = 0

		# if we're on the forward strand, use the "start" coordinate in forward direction
		if self.strand == 1 :

			leftBound = self.start - TF_RANGE_BEFORE
			rightBound = self.start + TF_RANGE_AFTER

		# if we're on the reverse strand, use the "end" coordinate in reverse direction
		elif self.strand == -1 :

			leftBound = self.end - TF_RANGE_AFTER
			rightBound = self.end + TF_RANGE_BEFORE

		# strand wasn't recognized
		else :
			print("ERROR: strand not recognized: \"%i\"" % self.strand)
			sys.exit(-1) 


		# sanity check on bounds
		if leftBound >= rightBound or max(leftBound, rightBound) <= 0 :
			print("ERROR: right/left bound exception: (%i,%i)" % (leftBound, rightBound))
			sys.exit(-1)

		# check if the TF is within the bounds :
		if TFHit.mid >= leftBound and TFHit.mid <= rightBound :
			return True

		else : return False

	def getRelativeHitPos(self, TFHit) :

		# if we're on the forward strand, use normal l->r direction from start
		if self.strand == 1 :

			return TFHit.mid - self.start

		# if we're on the reverse strand, use r->l direction from end
		elif self.strand == -1 :

			return self.end - TFHit.mid

	def getAbsoluteHitPos(self, TFHit) :

		# if we're on the forward strand, use normal l->r direction from start
		if self.strand == 1 :

			return np.abs(TFHit.mid - self.start)

		# if we're on the reverse strand, use r->l direction from end
		elif self.strand == -1 :

			return np.abs(TFHit.mid - self.end)

# class to hold a collection of upregulated genes
class GeneCollection :

	# read in the gene file, and store a list of the genes
	def __init__(self, pathToGeneFile) :

		"""
		expected schema for gene file (tsv) :
			Gene_WBID	Gene_Public_Name	Chr_Name	Strand	Start(bp)	End(bp)
		"""

		self.genes = []
		firstLine = False # keep track of whether we've parsed the header
		for line in open(pathToGeneFile, "rU") :

			if not firstLine and line[0] == "G" :
				firstLine = True
				continue
			try :
				l = line.strip().split("\t")
				WBID = l[0]
				PubName = l[1]
				strand = int(l[3])
				stPos = int(l[4])
				endPos = int(l[5])
				chrm = l[2]

				self.genes.append(Gene(chrm, strand, stPos, endPos, WBID, PubName))

			except :
				continue

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