groupCount = 1

with open("paralogGroups", "r") as f :
	for line in f :
		NEPGCount = 0
		l = line.strip().split()
		for geneID in l :
			if geneID.strip().split('_')[0] == "NEPG" :
				NEPGCount += 1
		print("group: %i #NEPG: %i total: %i" % (groupCount, NEPGCount, len(l)))
		groupCount += 1c