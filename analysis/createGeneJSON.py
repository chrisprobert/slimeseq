import json
from Bio import SeqIO, SearchIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML, Record
import os


links = []
linksD = {}
nodes = []
nodesByID = {}


with open("../inputSeqs/collections/4Nematocida/totalProteins.fa", "r") as prots :
	for record in SeqIO.parse(prots, "fasta") :

		thisIndex = len(nodes)
		thisID = record.id

		nodes.append({"name": thisID, "groups": [], "id": thisID})
		nodesByID[thisID] = thisIndex


groupCount = 0
with open("../orthoMCL_clusters/4Nematocida/paralogGroups", "r") as fin:
	for line in fin :
		groupCount += 1
		group = line.strip().split()

		for p in group :
			nodeIndex = nodesByID[p]
			nodes[nodeIndex]["groups"].append(groupCount)

myResults = os.listdir("../blastSearches/blastResults/")
for result in myResults :
	try :
		#with open("../blastResults/" + result, "r") as blastResultsHandle :
		#blastParser = NCBIXML.parse(blastResultsHandle)
		blastParser = SearchIO.parse("../blastSearches/blastResults/" + result, "blast-xml")
		for blastRecord in blastParser :
			queryIndex = nodesByID[blastRecord.id]
			for hit in blastRecord.hits :
				hit_id = hit[0].hit_id
				score = hit[0].bitscore_raw
				hitIndex = nodesByID[hit_id]
				id1 = "%s_%s" % (queryIndex, hitIndex)
				id2 = "%s_%s" % (hitIndex, queryIndex)
				if (id1 not in linksD) and (id2 not in linksD) :
					linksD[id1] = True
					linksD[id2] = True
					links.append({"source": queryIndex, "target": hitIndex, "value": score})
	except :
		print("Error on %s" % result)
		continue

with open("../blastResults/nodes.json", "w") as f :
	f.write(json.dumps(nodes))

with open("../blastResults/links.json", "w") as f :
	f.write(json.dumps(links))