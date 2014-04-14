import sys
import os


targetPOutputFile = "targetP_out.txt"

outputFile = "targetPSignals.txt"

with open(targetPOutputFile, "rU") as fin :
	with open(outputFile, "w") as fout :
		for line in f :
			lineSplit = line.strip().split()
			if len(lineSplit) == 7 :
				name = lineSplit[0]
				SP = float(lineSplit[3])
				Loc = lineSplit[5]
				RC = int(lineSplit[6])
				
				if Loc == "S" :
					fout.write(name + "\t" + str(SP) + "\t" + Loc + "\t" + RC + "\n")
				
				
