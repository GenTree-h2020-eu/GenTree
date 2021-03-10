#finalize HDplot results and create diagnostic plots
#christian rellstab 2020
#example on FS samples

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("geneplotter")
library(geneplotter)
library(data.table)
rm(list=ls(all=TRUE))

#change settings here
species <- "FS"
missing.results <- "Fsylvatica_GM_Oulu_filtered_v2._RENAMED.vcf.miss.lmiss"
HDplot.results <- "Fsylvatica_GM_Oulu_filtered_v2._RENAMED.depths"

#import SNP list with no. of missing data
missing <- read.table(missing.results, sep="\t", header=TRUE)
#missing data is given in haplotypes, has to becorrected
MISS_IND <- missing$N_MISS/2
#add SNP identifier
SNP.missing <- paste(missing$CHR,"_",missing$POS, sep="")
#combine the two vectors
missing2 <- cbind(SNP.missing, MISS_IND)

#import HD plot results (normally less SNPs than above)
HDplot <- read.table(HDplot.results, sep="\t", header=TRUE)
#add SNP identifier
SNP.HDplot <- paste(HDplot$contig,"_",HDplot$pos, sep="")
#create new table and remove unnecessary columns
HDplot2 <- cbind(SNP.HDplot, HDplot)
HDplot3 <- subset(HDplot2, select = -c(X, locus_ID))

#merge the two tables
table <- merge(missing2, HDplot3, by.x="SNP.missing", by.y="SNP.HDplot", all.x=T)

#calculate the correct H and AB
H <- as.numeric(table$num_hets)/(as.numeric(table$num_samples)-as.numeric(as.character(table$MISS_IND)))
AB <- pmax(table$depth_a,table$depth_b)/pmin(table$depth_a,table$depth_b)

#create and save final table
table2<-cbind(table,H,AB)
setnames(table2, old=c("ratio", "z"), new=c("RAF", "D"))
write.table(table2,paste(species, "_HDplot_results.txt", sep=""),sep="\t",quote=F, row.names=F)

#create a summary plot
#change thresholds here (this is just for plotting)
Hmax <- 0.6
RAFmin <- 0.2
RAFmax <- 0.8
Dmin <- -20
Dmax <- 20

#plot
par(mfrow=c(2,5))
hist(table2$H, breaks=50, xlab="H", main=species, cex.lab=1.2, cex.axis=1.2)
abline(v=Hmax, lty=2)
hist(table2$RAF, xlab="RAF [-]", breaks=50, main=species, cex.lab=1.2, cex.axis=1.2)
abline(v=RAFmax, lty=2)
abline(v=RAFmin, lty=2)
hist(table2$D, xlab="D", breaks=100, xlim=(c(-100,100)), main=species, cex.lab=1.2, cex.axis=1.2)
abline(v=Dmax, lty=2)
abline(v=Dmin, lty=2)
plot(table2$total_depth, table2$RAF, main = species, cex=0.5, pch=".", ylab="RAF", xlab="Depth")
abline(h=RAFmax, lty=2)
abline(h=RAFmin, lty=2)
plot(table2$total_depth, table2$D, main = species, cex=0.5, pch=".", ylab="D", xlab="Depth")
abline(h=Dmax, lty=2)
abline(h=Dmin, lty=2)
plot(table2$RAF, table2$D, main = species, cex=0.5, pch=".", ylab="D", xlab="RAF")
abline(v=RAFmax, lty=2)
abline(v=RAFmin, lty=2)
abline(h=Dmax, lty=2)
abline(h=Dmin, lty=2)
plot(table2$H, table2$RAF, main = species, cex=0.5, pch=".", ylab="RAF", xlab="H")
abline(h=RAFmax, lty=2)
abline(h=RAFmin, lty=2)
abline(v=Hmax, lty=2)
smoothScatter(table2$H, table2$RAF, xlab="H", ylab="RAF", nrpoints=0, main=species, cex.lab=1.2, cex.axis=1.2)
abline(h=RAFmax, lty=2)
abline(h=RAFmin, lty=2)
abline(v=Hmax, lty=2)
plot(table2$H, table2$D, main = species, ylim=(c(-100,100)), cex=0.5, pch=".", ylab="D", xlab="H")
abline(h=Dmax, lty=2)
abline(h=Dmin, lty=2)
abline(v=Hmax, lty=2)
smoothScatter(table2$H, table2$D, xlab="H", ylab="D", nrpoints=0, main=species, ylim=(c(-100,100)), cex.lab=1.2, cex.axis=1.2)
abline(h=RAFmax, lty=2)
abline(h=Dmax, lty=2)
abline(h=Dmin, lty=2)

#save as 5x10 landscape

#dev.off()
