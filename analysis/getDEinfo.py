import glob


PATH_TO_DIFFEXPR_FILES = "../input/diffexpr/"
DE_FILES = glob.glob(PATH_TO_DIFFEXPR_FILES + "*.results.txt")


def getAllUpregGeneWBIDs() :

	upregWBIDs = []

	for DEFileName in DE_FILES :
			for line in open( DEFileName ) :

				if line[0] == "l" : continue # skip header line

				l = line.strip().split()
				WBID = l[0].split(":")[0]
				logFC = float(l[2])
				PValue = float(l[3])
				FDR = float(l[4])

				# only taking upregulated genes:
				if logFC < 0 :

					if WBID not in upregWBIDs : upregWBIDs.append(WBID)

	return upregWBIDs


def getAllDownregGeneWBIDs() :

	downregWBIDs = []

	for DEFileName in DE_FILES :
			for line in open( DEFileName ) :

				if line[0] == "l" : continue # skip header line

				l = line.strip().split()
				WBID = l[0].split(":")[0]
				logFC = float(l[2])
				PValue = float(l[3])
				FDR = float(l[4])

				# only taking upregulated genes:
				if logFC > 0 :

					if WBID not in downregWBIDs : downregWBIDs.append(WBID)

	return downregWBIDs

def getUpDownIntersection() :
	inter = []
	up = getAllUpregGeneWBIDs()
	down = getAllDownregGeneWBIDs()
	for _id in up :
		if _id in down :
			inter.append(_id)
	return inter