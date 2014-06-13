#
#
#	file: getCorrelatedGeneList.py
#	date: 05/20/14
#	author: Chris Probert
#	
#	Purpose:
#	Generate a list of genes with nearby TF binding sites for a given TF.
#
#
#

import os, glob, sys
import numpy as np
import scipy.misc as scm
from dataModel import GFFfile, TFhit, Gene, GeneCollection, probabilityModel
import analysisSettings as settings

# path to the input data. All ".gff3" files are used.
PATH_TO_GFFS = "../input/modEncode/refinedPeaks/"
GFF_FILES = glob.glob(PATH_TO_GFFS + "*.gff3") + glob.glob(PATH_TO_GFFS + "*.GFF3")

# path to wormbase genome file
PATH_TO_GENEFILE = "../genome/wb_220bugfix_allGenes_stranded.txt"

# path to the list of genes from Grove paper
PATH_TO_GROVE_GENELIST = "../input/GroveMDL1Genes.txt"

# path to reference genome. Used for getting chromosome sizes.
PATH_TO_GENOME = "../genome/c_elegans.WS220.genomic.fa"

# path to directory containing differential expression tables
PATH_TO_DIFFEXPR_FILES = "../input/diffexpr/"
DE_FILES = glob.glob(PATH_TO_DIFFEXPR_FILES + "*.results.txt")

# window size for building null model of TF binding activity
WINDOW_SIZE = 200000

debug = True


def main() :

	print("loading analysis settings")
	TF_RANGE_BEFORE = settings.TF_RANGE_BEFORE
	TF_RANGE_AFTER = settings.TF_RANGE_AFTER
	stub = settings.stub

	print("reading gene metadata..")
	genes = GeneCollection(PATH_TO_GENEFILE)
	if debug: print("Num genes total: %i" % len(genes.genes))

	print("reading gene DE data..")
	genes = updateGenesDEData(genes, DE_FILES)
	if debug: print("Num DE genes total: %i" % len([g for g in genes.genes if g.DE == True]))

	print("reading gene Grove paper data..")
	genes = addGeneGroveOverlap(genes, PATH_TO_GROVE_GENELIST)
	if debug: print("Num Grove overlap genes total: %i" % len([g for g in genes.genes if g.GroveOverlap]))

	print("creating background probability model..")
	GFFs = [GFFfile(g) for g in GFF_FILES if "MDL-1" in g or "MDL1" in g ]
	allTFHits = []
	for GFF in GFFs :
		print("reading file %s" % GFF.filePath)
		for hit in GFF.hits : allTFHits.append(hit)

	with open(settings.Tables_Path + "DE_MDL1_%s_Intersect.csv" % stub, "w") as writer :

		writer.write("%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % \
			("Gene_Public_Name", \
			 "Gene_WB220ID", \
			 "DE_Pval", \
			 "DE_FDR", \
			 "DE_logFC", \
			 "TFhit_pval", \
			 "TFHit_MinAbsDistance", \
			 "TFHit_MinRelDistance", \
			 "TotalHits_in_range"))
		


		DEGenes = [g for g in genes.genes if g.DE == True]
		for myGene in DEGenes :

			# some quick sanity checks on this gene object
			# will make sure we don't join empty genes with empty TF sites
			if myGene.start == myGene.end : continue
			if min(myGene.start, myGene.end) < 1 : continue
			if len(myGene.chrm.strip()) < 1 : continue

			myGeneHits = [ hit for hit in allTFHits if \
							myGene.checkTFHitValid(hit, TF_RANGE_BEFORE, TF_RANGE_AFTER) ]
			
			if len(myGeneHits) >= 1 :

				minAbsDistance = np.min([ myGene.getAbsoluteHitPos(hit) for hit in myGeneHits ])
				minRelDistance = [ myGene.getRelativeHitPos(hit) for hit in myGeneHits if \
											myGene.getAbsoluteHitPos(hit) == minAbsDistance ][0]

				writer.write("%s,%s,%.8f,%.8f,%.8f,%.8f,%i,%i,%i\n" % \
					(myGene.publicName, \
					 myGene.WBGeneID, \
					 myGene.PValue, \
					 myGene.FDR,\
					 myGene.logFC,\
					 myGeneHits[0].pval,\
					 minAbsDistance,\
					 minRelDistance,\
					 len(myGeneHits)))

	with open(settings.Tables_Path + "DE_Grove-MDL1_Intersect.csv", "w") as writer :
		
		writer.write("%s,%s,%s,%s,%s\n" % \
			("Gene_Public_Name", \
			 "Gene_WB220ID", \
			 "DE_Pval", \
			 "DE_FDR", \
			 "DE_logFC" ))
		
		for myGene in [ g for g in genes.genes if g.GroveOverlap and g.DE ] :

			writer.write("%s,%s,%.8f,%.8f,%.8f\n" % \
				(myGene.publicName, \
				 myGene.WBGeneID, \
				 myGene.PValue, \
				 myGene.FDR,\
				 myGene.logFC ))

	with open(settings.Tables_Path + "DE_MDL1_%s_Grove-MDL1_3-way-Intersect.csv" % stub, "w") as writer :

		writer.write("%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % \
			("Gene_Public_Name", \
			 "Gene_WB220ID", \
			 "DE_Pval", \
			 "DE_FDR", \
			 "DE_logFC", \
			 "TFhit_pval", \
			 "TFHit_MinAbsDistance", \
			 "TFHit_MinRelDistance", \
			 "TotalHits_in_range"))
		


		DEGenes = [g for g in genes.genes if g.GroveOverlap and g.DE == True]
		for myGene in DEGenes :

			# some quick sanity checks on this gene object
			# will make sure we don't join empty genes with empty TF sites
			if myGene.start == myGene.end : continue
			if min(myGene.start, myGene.end) < 1 : continue
			if len(myGene.chrm.strip()) < 1 : continue

			myGeneHits = [ hit for hit in allTFHits if \
							myGene.checkTFHitValid(hit, TF_RANGE_BEFORE, TF_RANGE_AFTER) ]
			
			if len(myGeneHits) >= 1 :

				minAbsDistance = np.min([ myGene.getAbsoluteHitPos(hit) for hit in myGeneHits ])
				minRelDistance = [ myGene.getRelativeHitPos(hit) for hit in myGeneHits if \
											myGene.getAbsoluteHitPos(hit) == minAbsDistance ][0]

				writer.write("%s,%s,%.8f,%.8f,%.8f,%.8f,%i,%i,%i\n" % \
					(myGene.publicName, \
					 myGene.WBGeneID, \
					 myGene.PValue, \
					 myGene.FDR,\
					 myGene.logFC,\
					 myGeneHits[0].pval,\
					 minAbsDistance,\
					 minRelDistance,\
					 len(myGeneHits)))

	with open(settings.Tables_Path + "AllGenes_MDL1_%s_Intersect.csv" % stub, "w") as writer :
		writer.write("%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % \
			("Gene_Public_Name",\
			 "Gene_WB220ID",\
			 "DE_Pval",\
			 "DE_FDR",\
			 "DE_logFC",\
			 "TFhit_pval",\
			 "TFHit_MinAbsDistance",\
			 "TFHit_MinRelDistance",\
			 "TotalHits_in_range"))

		DEGenes = [g for g in genes.genes]

		for myGene in DEGenes :

			# some quick sanity checks on this gene object
			# will make sure we don't join empty genes with empty TF sites
			if myGene.start == myGene.end : continue
			if min(myGene.start, myGene.end) < 1 : continue
			if len(myGene.chrm.strip()) < 1 : continue

			myGeneHits = [ hit for hit in allTFHits if \
								myGene.checkTFHitValid(hit, TF_RANGE_BEFORE, TF_RANGE_AFTER) ]
			
			if len(myGeneHits) >= 1 :

				minAbsDistance = np.min([ np.abs(hit.mid - myGene.start) for hit in myGeneHits ])
				minRelDistance = [ (hit.mid - myGene.start) for hit in myGeneHits if \
											np.abs(hit.mid - myGene.start) == minAbsDistance ][0]

				writer.write("%s,%s,%.8f,%.8f,%.8f,%.8f,%i,%i,%i\n" % \
					(myGene.publicName,\
					 myGene.WBGeneID,\
					 myGene.PValue,\
					 myGene.FDR,\
					 myGene.logFC,\
					 myGeneHits[0].pval,\
					 minAbsDistance,\
					 minRelDistance,\
					 len(myGeneHits)))

	with open(settings.Tables_Path + "AllGenes_Grove-MDL1_Intersect.csv", "w") as writer :

		writer.write("%s,%s,%s,%s,%s\n" % \
			("Gene_Public_Name", \
			 "Gene_WB220ID", \
			 "DE_Pval", \
			 "DE_FDR", \
			 "DE_logFC" ))
		
		for myGene in [ g for g in genes.genes if g.GroveOverlap ] :

			writer.write("%s,%s,%.8f,%.8f,%.8f\n" % \
				(myGene.publicName, \
				 myGene.WBGeneID, \
				 myGene.PValue, \
				 myGene.FDR,\
				 myGene.logFC ))

	with open(settings.Tables_Path + "AllGenes_MDL1_%s_Grove-MDL1_3-way-Intersect.csv" % stub, "w") as writer :
		writer.write("%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % \
			("Gene_Public_Name",\
			 "Gene_WB220ID",\
			 "DE_Pval",\
			 "DE_FDR",\
			 "DE_logFC",\
			 "TFhit_pval",\
			 "TFHit_MinAbsDistance",\
			 "TFHit_MinRelDistance",\
			 "TotalHits_in_range"))

		DEGenes = [g for g in genes.genes if g.GroveOverlap ]

		for myGene in DEGenes :

			# some quick sanity checks on this gene object
			# will make sure we don't join empty genes with empty TF sites
			if myGene.start == myGene.end : continue
			if min(myGene.start, myGene.end) < 1 : continue
			if len(myGene.chrm.strip()) < 1 : continue

			myGeneHits = [ hit for hit in allTFHits if \
								myGene.checkTFHitValid(hit, TF_RANGE_BEFORE, TF_RANGE_AFTER) ]
			
			if len(myGeneHits) >= 1 :

				minAbsDistance = np.min([ np.abs(hit.mid - myGene.start) for hit in myGeneHits ])
				minRelDistance = [ (hit.mid - myGene.start) for hit in myGeneHits if \
											np.abs(hit.mid - myGene.start) == minAbsDistance ][0]

				writer.write("%s,%s,%.8f,%.8f,%.8f,%.8f,%i,%i,%i\n" % \
					(myGene.publicName,\
					 myGene.WBGeneID,\
					 myGene.PValue,\
					 myGene.FDR,\
					 myGene.logFC,\
					 myGeneHits[0].pval,\
					 minAbsDistance,\
					 minRelDistance,\
					 len(myGeneHits)))


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
					if genes.genes[i].DE == False or (genes.genes[i].FDR > FDR) :
						genes.genes[i].addDEInfo(logFC, PValue, FDR)

	return genes

def addGeneGroveOverlap(genes, groveFilePath) :

	groveGenes = [ g.strip() for g in open(groveFilePath) ]

	for g in genes.genes :
		if g.WBGeneID in groveGenes : g.addGroveOverlap()

	return genes


if __name__ == "__main__" : main()