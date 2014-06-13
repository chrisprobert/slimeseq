import matplotlib.pyplot as plt
from getDEinfo import getAllUpregGeneWBIDs, getAllDownregGeneWBIDs
import analysisSettings as settings
import glob


stub = settings.stub
TFRangeBefore = settings.TF_RANGE_BEFORE
TFRangeAfter = settings.TF_RANGE_AFTER
# historgram of MDL-1 binding sites relative to all gene TSSs

distances = [ int( line.strip().split(',')[7] ) for line in \
				[ l for l in open(settings.Tables_Path + "AllGenes_MDL1_%s_Intersect.csv"%stub) ][1:] ]
plt.clf()
plt.hist(distances, bins=range(-TFRangeBefore,TFRangeAfter+50,50))
plt.axvline(x=-700, color='r')
plt.axvline(x=100, color='r')

plt.title("Position of nearest MDL-1 binding site - all genes")
plt.xlabel("Relative position from TSS")
plt.ylabel("Number of genes (bin width = 50)")
plt.savefig(settings.Figures_Path + "TSS_relative_dists_AllGenes_%s"%stub)



# historgram of MDL-1 binding sites relative to DE gene TSSs

distances = [ int( line.strip().split(',')[7] ) for line in \
				[ l for l in open(settings.Tables_Path + "DE_MDL1_%s_Intersect.csv"%stub) ][1:] ]
plt.clf()
plt.hist(distances, bins=range(-TFRangeBefore,TFRangeAfter+100,100))

plt.title("Position of nearest MDL-1 binding site - DE genes only")
plt.xlabel("Relative position from TSS")
plt.ylabel("Number of genes (bin width = 100)")
plt.savefig(settings.Figures_Path + "TSS_relative_dists_DE_%s"%stub)



# historgram of MDL-1 binding sites relative to upreg-DE gene TSSs

upregWBIDS = getAllUpregGeneWBIDs()
window_WBIDS = [ line.strip().split(',')[1] for line in \
					[ l for l in open( settings.Tables_Path + "DE_MDL1_%s_Intersect.csv"%stub )][1:] ]
with open(settings.Tables_Path + "WBIDs_upreg_%s_intersect.txt"%stub, "w") as fout :
	for WBID in upregWBIDS :
		if WBID in window_WBIDS : fout.write("%s\n" % WBID)

distances = [ int( line.strip().split(',')[7] ) for line in \
				[ l for l in open(settings.Tables_Path + "DE_MDL1_%s_Intersect.csv"%stub) ][1:] \
				if line.strip().split(',')[1] in upregWBIDS ]

plt.clf()
plt.hist(distances, bins=range(-TFRangeBefore,TFRangeAfter+150,150))

plt.title("Position of nearest MDL-1 binding site - upregulated genes only")
plt.xlabel("Relative position from TSS")
plt.ylabel("Number of genes (bin width = 150)")
plt.savefig(settings.Figures_Path + "TSS_relative_dists_upreg_%s"%stub)



# historgram of MDL-1 binding sites relative to downreg-DE gene TSSs

downregWBIDS = getAllDownregGeneWBIDs()
window_WBIDS = [ line.strip().split(',')[1] for line in \
					[ l for l in open( settings.Tables_Path + "DE_MDL1_%s_Intersect.csv"%stub )][1:] ]
with open(settings.Tables_Path + "WBIDs_downreg_%s_intersect.txt"%stub, "w") as fout :
	for WBID in downregWBIDS :
		if WBID in window_WBIDS : fout.write("%s\n" % WBID)

distances = [ int( line.strip().split(',')[7] ) for line in \
				[ l for l in open(settings.Tables_Path + "DE_MDL1_%s_Intersect.csv"%stub) ][1:] \
				if line.strip().split(',')[1] in downregWBIDS ]
plt.clf()
plt.hist(distances, bins=range(-TFRangeBefore,TFRangeAfter+150,150))

plt.title("Position of nearest MDL-1 binding site - downregulated genes only")
plt.xlabel("Relative position from TSS")
plt.ylabel("Number of genes (bin width = 150)")
plt.savefig(settings.Figures_Path + "TSS_relative_dists_downreg_%s"%stub)