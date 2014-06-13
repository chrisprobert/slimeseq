
# Paths to directories where tables and figures are stored

Tables_Path = "./tables/"
Figures_Path = "./figures/"

# binding range for how far a TF hit can be from a gene (in bp)
TF_RANGE_BEFORE = 700
TF_RANGE_AFTER = 100

# stub to identify files with given upstream/downstream ranges
stub = "-%i_%i" % (TF_RANGE_BEFORE, TF_RANGE_AFTER)