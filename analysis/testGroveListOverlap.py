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
from getDEinfo import getAllUpregGeneWBIDs, getAllDownregGeneWBIDs
import analysisSettings as settings

# path to the input data. All ".gff3" files are used.
PATH_TO_GFFS = "../input/modEncode/refinedPeaks/"
GFF_FILES = glob.glob(PATH_TO_GFFS + "*.gff3") + glob.glob(PATH_TO_GFFS + "*.GFF3")

# path to wormbase genome file
PATH_TO_GENEFILE = "../genome/wb_220bugfix_allGenes_stranded.txt"

# path to the list of genes from Grove paper
PATH_TO_GENELIST = "../input/Grove2009Cell_sup7.csv"

# path to directory containing differential expression tables
PATH_TO_DIFFEXPR_FILES = "../input/diffexpr/"
DE_FILES = glob.glob(PATH_TO_DIFFEXPR_FILES + "*.results.txt")

# window size for building null model of TF binding activity
WINDOW_SIZE = 200000

debug = True


def main() :

	print("loading analysis settings")
	stub = settings.stub

	print("reading gene metadata..")
	genes = GeneCollection(PATH_TO_GENEFILE)
	if debug: print("Num genes total: %i" % len(genes.genes))

	print("reading gene DE data..")
	genes = updateGenesDEData(genes, DE_FILES)
	if debug: print("Num DE genes total: %i" % len([g for g in genes.genes if g.DE == True]))

	print("reading GFF files..")
	GFFs = [GFFfile(g) for g in GFF_FILES if "MDL-1" in g or "MDL1" in g ]
	allTFHits = []
	for GFF in GFFs :
		print("reading file %s" % GFF.filePath)
		for hit in GFF.hits : allTFHits.append(hit)

	print("reading grove gene list..")
	# we know 16 is the column index for MDL-1 (checked)
	# the first row is header (label)
	GroveMDL1Genes = [ l.split(',')[16] for l in open("../input/Grove2009Cell_sup7.csv") ]
	if debug: print("Grove table column header: %s" % GroveMDL1Genes[-1])
	# now remove the column header
	GroveMDL1Genes = GroveMDL1Genes[1:]
	if debug: print("Num Grove genes total: %i" % len(GroveMDL1Genes))



	# write a list of all Grove paper genes. Report an error if gene not in our list.
	with open(settings.Tables_Path + "GrovePaper-MDL1_AllGenes.csv" , "w") as writer :

		writer.write("%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % \
			("Gene_Public_Name", \
			 "Gene_WB220ID", \
			 "DE_Pval", \
			 "DE_FDR", \
			 "DE_logFC" ))
		
		for WBGeneID in GroveMDL1Genes :
			
			myGene = [ g for g in genes.genes if g.WBGeneID == WBGeneID ][0]

			writer.write("%s,%s,%.8f,%.8f,%.8f\n" % \
				(myGene.publicName, \
				 myGene.WBGeneID, \
				 myGene.PValue, \
				 myGene.FDR,\
				 myGene.logFC ))

	# write a list of Grove paper genes that are DE upon infection.
	with open(settings.Tables_Path + "GrovePaper-MDL1_DEGenes_Intersect.csv" % stub, "w") as writer :
		writer.write("%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % \
			("Gene_Public_Name",\
			 "Gene_WB220ID",\
			 "DE_Pval",\
			 "DE_FDR",\
			 "DE_logFC" ))

		DEGenes = [g for g in genes.genes if g.DE == True]

		for myGene in DEGenes :
			if myGene.WBGeneID in GroveMDL1Genes :

				writer.write("%s,%s,%.8f,%.8f,%.8f,%.8f,%i,%i,%i\n" % \
					(myGene.publicName,\
					 myGene.WBGeneID,\
					 myGene.PValue,\
					 myGene.FDR,\
					 myGene.logFC ))

	# write a list of Grove paper genes that are upregulated upon infection.
	with open(settings.Tables_Path + "GrovePaper-MDL1_Upregulated_Intersect.csv" % stub, "w") as writer :
		writer.write("%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % \
			("Gene_Public_Name",\
			 "Gene_WB220ID",\
			 "DE_Pval",\
			 "DE_FDR",\
			 "DE_logFC" ))

		targetGenes = getAllUpregGeneWBIDs()
		DEGenes = [ g for g in genes.genes if g.WBGeneID in targetGenes ]

		for myGene in DEGenes :
			if myGene.WBGeneID in GroveMDL1Genes :

				writer.write("%s,%s,%.8f,%.8f,%.8f,%.8f,%i,%i,%i\n" % \
					(myGene.publicName,\
					 myGene.WBGeneID,\
					 myGene.PValue,\
					 myGene.FDR,\
					 myGene.logFC ))

	# write a list of Grove paper genes that are upregulated upon infection.
	with open(settings.Tables_Path + "GrovePaper-MDL1_Downregulated_Intersect.csv" % stub, "w") as writer :
		writer.write("%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % \
			("Gene_Public_Name",\
			 "Gene_WB220ID",\
			 "DE_Pval",\
			 "DE_FDR",\
			 "DE_logFC" ))

		targetGenes = getAllDownregGeneWBIDs()
		DEGenes = [ g for g in genes.genes if g.WBGeneID in targetGenes ]

		for myGene in DEGenes :
			if myGene.WBGeneID in GroveMDL1Genes :

				writer.write("%s,%s,%.8f,%.8f,%.8f,%.8f,%i,%i,%i\n" % \
					(myGene.publicName,\
					 myGene.WBGeneID,\
					 myGene.PValue,\
					 myGene.FDR,\
					 myGene.logFC ))




if __name__ == "__main__" : main()