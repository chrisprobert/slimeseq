import sys
import os


targetPOutputFile = "signalPSignals.txt"

outputFile = "totalProteins_Signalp_out"

with open(outputFile, "rU") as fin :
	with open(targetPOutputFile, "w") as fout :
		for line in fin :
			lineSplit = line.strip().split()
			if len(lineSplit) == 12 :
				if lineSplit[9] == "Y" :
					name = lineSplit[0]
					network = lineSplit[11]
					Dmax = lineSplit[8]
					fout.write("\n" + name + "\t" + network + "\t" + Dmax)
				
				
