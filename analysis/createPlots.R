setwd("/Users/biotools/MicrosporidiaPhyloComp/NEPG_GroupAnalysis/output_data/")

require("MASS")
require("ggplot2")
require("limma")



clusterSizes <- read.table(file="resultsByCluster.csv", header=TRUE, sep=",")
clusters2NEPG <- clusterSizes[clusterSizes$NEPG > 2,]
clusters2TPHits <- clusterSizes[clusterSizes$targetPHits > 2,]
clusters2NEPGTPHits <- clusters2NEPG[clusters2NEPG$targetPHits > 1,]
clustersHEPOverlap <- clusterSizes[clusterSizes$proteomicsOverlap > 0,]
clustersHMMhits <- clusterSizes[clusterSizes$HMMerHits > 0,]


NEPGclusters <- read.table(file="resultsByCluster_NEPGOnly.csv", header=TRUE, sep=",")
NEPGclusters10NEPG <- NEPGclusters[NEPGclusters$NEPG > 2,]
NEPGclustersTPHits <- NEPGclusters[NEPGclusters$targetPHits > 1,]
NEPGclustersHMMhits <- NEPGclusters[NEPGclusters$HMMerHits > 1,]


NGF1 <- read.table(file="../NGF1.txt")
NGF2 <- read.table(file="../NGF2.txt")
NGF3 <- read.table(file="../NGF3.txt")
prots <- read.table(file="../HostExposedProts")
NEPG_Genome_Size <- 2661

print("Size of NGF1:")
print(nrow(NGF1))
print("Size of NGF2:")
print(nrow(NGF2))
print("Size of NGF3:")
print(nrow(NGF3))
print("Size of Host-Exposed Protein List:")
print(nrow(prots))


print("Clusters with NEPG > 2:")
print(clusters2NEPG)

print("Clusters with NEPG > 2 and TargetP hits > 1:")
print(clusters2NEPGTPHits)

print("Clusters with proteomics overlap:")
print(clustersHEPOverlap)

print("NEPG-only clusters with > 10 proteins:")
print(NEPGclusters10NEPG)

print("NEPG-only clusters with TargetP hits > 1:")
print(NEPGclustersTPHits)


#
# Scatterplots of targetP hits vs. cluster size
#
ggplot(NEPGclustersTPHits, aes(NEPG,targetPHits)) + geom_point() + geom_text(aes(label=clusterNum), hjust= -1) + xlab("Num. NEPG Proteins") + ylab("Num. TargetP Hits") + ggtitle("TargetP Hits by NEPG-only Cluster Size") + ggsave("Scatter_TargetPHitsNEPG_NEPGOnly.pdf")

ggplot(clusters2NEPG, aes(NEPG,targetPHits)) + geom_point() + geom_text(aes(label=clusterNum), hjust= -1) + xlab("Num. NEPG Proteins") + ylab("Num. TargetP Hits") + ggtitle("TargetP Hits by Cluster Size in Clusters with > 2 NEPG Proteins") + ggsave("Scatter_TargetPHitsNEPG_AllClusters.pdf")


#
# Barcharts of Hmmer hits / targetP hits by NEPG-only cluster 
#
# NEPG-only clusters
NEPGclustersHMMhits$clusterNum <- as.character(NEPGclustersHMMhits$clusterNum)
NEPGclustersHMMhits$clusterNum <- factor(NEPGclustersHMMhits$clusterNum, levels=unique(NEPGclustersHMMhits$clusterNum), ordered=TRUE)
NEPGclustersHMMhits$propHMMHitsTPHits <- NEPGclustersHMMhits$HMMtargetPHits / NEPGclustersHMMhits$HMMerHits
ggplot(NEPGclustersHMMhits, aes(x=clusterNum, y=HMMerHits, fill=propHMMHitsTPHits)) + geom_bar(stat="identity") + ggtitle("HMMer Hits by cluster in NEPG-only clusters with targetP hits") +  geom_text(aes(label=paste("TP:",HMMtargetPHits)), size=2.5, vjust=-.5) + xlab("NEPG-only cluster number") + ylab("num. NEPG-only HMMer hits") + ggsave("bar_HMMerTPHits_NEPGonlyClusters.pdf")

# all clusters
clustersHMMhits$clusterNum <- as.character(clustersHMMhits$clusterNum)
clustersHMMhits$clusterNum <- factor(clustersHMMhits$clusterNum, levels=unique(clustersHMMhits$clusterNum), ordered=TRUE)
clustersHMMhits$propHMMHitsTPHits <- clustersHMMhits$HMMtargetPHits / clustersHMMhits$HMMerHits
clustersHMMhits <- clustersHMMhits[clustersHMMhits$HMMerHits > 25,]
clustersHMMhits <- clustersHMMhits[clustersHMMhits$propHMMHitsTPHits > 0.25,]
ggplot(clustersHMMhits, aes(x=clusterNum, y=HMMerHits, fill=propHMMHitsTPHits)) + geom_bar(stat="identity") + ggtitle("HMMer Hits by cluster in clusters with targetP hits") +  geom_text(aes(label=paste("TP:",HMMtargetPHits)), size=2.5, vjust=-.5) + xlab("cluster number") + ylab("num. HMMer hits") + theme(text = element_text(size=8)) + ggsave("bar_HMMerTPHits_allClusters.pdf", width=12)


#
# bar charts showing number of targetP hits per cluster for NEPG-only clusters (also by P-value of TP hits ratio)
#
propTPHits <- NEPGclustersTPHits[NEPGclustersTPHits$total > 5,c("clusterNum", "total", "targetPHits")]; propTPHits$clusterNum <- paste("cluster", propTPHits$clusterNum, sep=" "); propTPHits <- rbind(propTPHits, c("genome", as.numeric(2661), as.numeric(566))); propTPHits <- transform(propTPHits, total = as.numeric(total)); propTPHits <- transform(propTPHits, targetPHits = as.numeric(targetPHits)); propTPHits$prop <- propTPHits$targetPHits / propTPHits$total

for(i in 1:nrow(propTPHits)) propTPHits$pval[i] <- (prop.test(x=c(propTPHits$targetPHits[i],566),n=c(propTPHits$total[i],2661)))$p.val;
for(i in 1:nrow(propTPHits)) propTPHits$NLpval[i] <- -log(propTPHits$pval[i]);

ggplot(propTPHits, aes(clusterNum, prop)) + geom_bar(stat="identity") + geom_text(aes(label=paste(paste("n=", targetPHits), "/", total)), vjust=-.5, size=4) + xlab("NEPG-only Protein Family Cluster") + ylab("Proportion of proteins with TargetP hits") + ggtitle("Proportion of NEPG proteins with TargetP hits by OrthoMCL cluster") + ggsave("bar_ProportionTargetPHitsNEPGonly.pdf")

ggplot(propTPHits, aes(clusterNum, NLpval)) + geom_bar(stat="identity") + geom_text(aes(label=paste("p=", sprintf('%3.1e',pval)), vjust=-.5, size=4)) + xlab("NEPG-only Protein Family Cluster") + ylab("-Log(Pval) for propotion of proteins with TargetP hits") + ggtitle("P-Values of TargetP hits by OrthoMCL cluster") + ggsave("bar_PValTargetPHitsNEPGonly.pdf")

#
# bar charts showing number of targetP hits per cluster for ALL clusters (also by P-value of TP hits ratio)
#
propTPHitsAll <- clusters2TPHits[clusters2TPHits$total > 5,c("clusterNum", "total", "targetPHits")];
propTPHitsAll <- propTPHitsAll[propTPHitsAll$targetPHits > 5,];  
propTPHitsAll$clusterNum <- paste("cluster", propTPHitsAll$clusterNum, sep=" "); propTPHitsAll <- rbind(propTPHitsAll, c("genome", as.numeric(2661), as.numeric(566))); propTPHitsAll <- transform(propTPHitsAll, total = as.numeric(total)); propTPHitsAll <- transform(propTPHitsAll, targetPHits = as.numeric(targetPHits)); propTPHitsAll$prop <- propTPHitsAll$targetPHits / propTPHitsAll$total

for(i in 1:nrow(propTPHitsAll)) propTPHitsAll$pval[i] <- (prop.test(x=c(propTPHitsAll$targetPHits[i],566),n=c(propTPHitsAll$total[i],2661)))$p.val;
for(i in 1:nrow(propTPHitsAll)) propTPHitsAll$NLpval[i] <- -log(propTPHitsAll$pval[i]);

ggplot(propTPHitsAll, aes(clusterNum, prop)) + geom_bar(stat="identity") + geom_text(aes(label=paste(paste("n=", targetPHits), "/", total)), vjust=-.5, size=4) + xlab("All Species Protein Family Cluster") + ylab("Proportion of proteins with TargetP hits") + ggtitle("Proportion of NEPG proteins with TargetP hits by OrthoMCL cluster") + ggsave("bar_ProportionTargetPHitsAllClusts.pdf", width=12)

ggplot(propTPHitsAll, aes(clusterNum, NLpval)) + geom_bar(stat="identity") + geom_text(aes(label=paste("p=", sprintf('%3.1e',pval)), vjust=-.5, size=4)) + xlab("All Species Protein Family Cluster") + ylab("-Log(Pval) for propotion of proteins with TargetP hits") + ggtitle("P-Values of TargetP hits by OrthoMCL cluster") + ggsave("bar_PValTargetPHitsAllClusts.pdf", width=12)


#
# bar charts showing overall cluster sizes colored by TP hits ratio
#
propTPHitsAll <- clusterSizes[clusterSizes$total > 5,c("clusterNum", "total", "targetPHits")]; 
propTPHitsAll <- transform(propTPHitsAll, total = as.numeric(total))
propTPHitsAll <- transform(propTPHitsAll, targetPHits = as.numeric(targetPHits))
propTPHitsAll$propTargetPHits <- propTPHitsAll$targetPHits / propTPHitsAll$total
ggplot(propTPHitsAll, aes(clusterNum, total, fill= propTargetPHits)) + geom_bar(stat="identity") + xlab("OrthoMCL Cluster") + ylab("Cluster size") + ggtitle("Nematocida-wide OrthoMCL cluster size and prop. targetP hits") + ggsave("bar_overallClustersTP_Size.pdf", width=12)



protStats <- read.table(file="proteinLevelStats.csv", sep=",", header=TRUE)
NEPG_Group_1 <- (protStats$group1 == 1)
NEPG_Group_3 <- (protStats$group3 == 1)
NEPG_Group_10 <- (protStats$group10 == 1)
NEPG_Group_22 <- (protStats$group22 == 1)
NGF1 <- (protStats$NGF1Overlap == 1)
NGF2 <- (protStats$NGF2Overlap == 1)
NGF3 <- (protStats$NGF3Overlap == 1)
Host_Exposed <- (protStats$proteomicsOverlap == 1)


pdf("VennDiag_NEPGGroup1_NGF1.pdf")
vennDiagram(vennCounts(cbind(NEPG_Group_1, NGF1)), main="Identity Overlap: NEPG_Group1 & NGF1")
dev.off()

pdf("VennDiag_NEPGGroup3_NGF2.pdf")
vennDiagram(vennCounts(cbind(NEPG_Group_3, NGF2)), main="Identity Overlap: NEPG_Group3 & NGF2")
dev.off()

pdf("VennDiag_NEPGGroup10_NGF3.pdf")
vennDiagram(vennCounts(cbind(NEPG_Group_10, NGF3)), main="Identity Overlap: NEPG_Group10 & NGF3")
dev.off()

pdf("VennDiag_NEPGGroup22_NGF3.pdf")
vennDiagram(vennCounts(cbind(NEPG_Group_22, NGF3)), main="Identity Overlap: NEPG_Group22 & NGF3")
dev.off()


ggplot(protStats, aes(factor(group1), length)) + geom_violin() + ggtitle("Distribution of protein lengths") + xlab("Frequency by OrthoMCL Group 1 membership") + ylab("Protein length (aa)") + scale_x_discrete(labels=c("genome","group 1")) + ggsave("Violin_LengthDistr_Group1.pdf")

ggplot(protStats, aes(factor(group2), length)) + geom_violin() + ggtitle("Distribution of protein lengths") + xlab("Frequency by OrthoMCL Group 2 membership") + ylab("Protein length (aa)") + scale_x_discrete(labels=c("genome","group 2")) + ggsave("Violin_LengthDistr_Group2.pdf")

ggplot(protStats, aes(factor(ProteomicsOverlap), length)) + geom_violin() + ggtitle("Distribution of protein lengths") + xlab("Frequency by Host Exposure") + ylab("Protein length (aa)") + scale_x_discrete(labels=c("genome","host exposed")) + ggsave("Violin_LengthDistr_HostExposure.pdf")



# make a scatterplot based on kmer frequency

km <- read.table(file="kmerMatrix_Group11.csv", header=TRUE, row.names=1, sep=",")
cleankm <- km[,!(names(km) %in% c("group"))]
cleankm <- sapply(cleankm, as.numeric); cleankmT <- t(cleankm)
d <- matrix(nrow=ncol(cleankmT), ncol=ncol(cleankmT))
for(i in 1:nrow(d))
	for(j in 1:nrow(d))
		d[i,j] <- 1 - cov(cleankmT[,i], cleankmT[,j]);
d[d==0] <- 0.00000001
fit <- as.data.frame((isoMDS(d, k=2))$points)
fit$labels <- km$group
ggplot(fit, aes(V1,V2)) + geom_point(aes(colour=labels)) + ggsave("scatter_kmer_PCA_group11.pdf")


km <- read.table(file="kmerMatrix_Group1.csv", header=TRUE, row.names=1, sep=",")
cleankm <- km[,!(names(km) %in% c("group"))]
cleankm <- sapply(cleankm, as.numeric); cleankmT <- t(cleankm)
d <- matrix(nrow=ncol(cleankmT), ncol=ncol(cleankmT))
for(i in 1:nrow(d))
	for(j in 1:nrow(d))
		d[i,j] <- 1 - cov(cleankmT[,i], cleankmT[,j]);
d[d==0] <- 0.00000001
fit <- as.data.frame((isoMDS(d, k=2))$points)
fit$labels <- km$group
ggplot(fit, aes(V1,V2)) + geom_point(aes(colour=labels)) + ggsave("scatter_kmer_PCA_group1.pdf")

km <- read.table(file="kmerMatrix_Group2.csv", header=TRUE, row.names=1, sep=",")
cleankm <- km[,!(names(km) %in% c("group"))]
cleankm <- sapply(cleankm, as.numeric); cleankmT <- t(cleankm)
d <- matrix(nrow=ncol(cleankmT), ncol=ncol(cleankmT))
for(i in 1:nrow(d))
	for(j in 1:nrow(d))
		d[i,j] <- 1 - cov(cleankmT[,i], cleankmT[,j]);
d[d==0] <- 0.00000001
fit <- as.data.frame((isoMDS(d, k=2))$points)
fit$labels <- km$group
ggplot(fit, aes(V1,V2)) + geom_point(aes(colour=labels)) + ggsave("scatter_kmer_PCA_group2.pdf")

km <- read.table(file="kmerMatrix_Group3.csv", header=TRUE, row.names=1, sep=",")
cleankm <- km[,!(names(km) %in% c("group"))]
cleankm <- sapply(cleankm, as.numeric); cleankmT <- t(cleankm)
d <- matrix(nrow=ncol(cleankmT), ncol=ncol(cleankmT))
for(i in 1:nrow(d))
	for(j in 1:nrow(d))
		d[i,j] <- 1 - cov(cleankmT[,i], cleankmT[,j]);
d[d==0] <- 0.00000001
fit <- as.data.frame((isoMDS(d, k=2))$points)
fit$labels <- km$group
ggplot(fit, aes(V1,V2)) + geom_point(aes(colour=labels)) + ggsave("scatter_kmer_PCA_group3.pdf")































































