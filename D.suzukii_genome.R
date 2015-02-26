contigStats <- function(N, reflength, style="ggplot2", pch=20, xlab="Percentage of Assembly Covered by Contigs of Size >=Y", ylab="Contig Size [bp]", main="Cumulative Length of Contigs", sizetitle=14, sizex=12, sizey=12, sizelegend=9, xlim, ylim) {
        ## Compute cumulative length vectors for contig sets
        Nl <- lapply(names(N), function(x) rev(sort(N[[x]]))); names(Nl) <- names(N)
        Nlcum <- lapply(names(Nl), function(x) cumsum(Nl[[x]])); names(Nlcum) <- names(Nl)
            
        ## Compute N50 values
        N50 <- sapply(seq(along=N), function(x) Nl[[x]][which(Nlcum[[x]] - reflength[x]/2 >= 0)[1]]); names(N50) <- names(N)
	
        ## Return only data (no plot)
        if(style=="data") {
                N75 <- sapply(seq(along=N), function(x) Nl[[x]][which(Nlcum[[x]] - reflength[x] * 0.75 >= 0)[1]]); names(N50) <- names(N)
                N25 <- sapply(seq(along=N), function(x) Nl[[x]][which(Nlcum[[x]] - reflength[x] * 0.25 >= 0)[1]]); names(N50) <- names(N)
                stats <- cbind(N25, N50, N75, Longest=sapply(N, max), Mean=sapply(N, mean), Median=sapply(N, median), Shortest=sapply(N, min), N_Contigs=sapply(N, length))
                return(c(Nlcum, Contig_Stats=list(stats)))
        }
        ## Plot cumulative contig length with base graphics
	if(style=="base") {
            if(missing(xlim)) xlim <- c(0, 100)
            if(missing(ylim)) ylim <- c(0, max(unlist(N))) 
            split.screen(c(1,1))
            for(i in seq(along=Nl)) {
                    if(i==1) {
                            plot(Nlcum[[i]]/reflength[[i]] * 100, Nl[[i]], col=i, pch=pch, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=main)  
                    } 
                    screen(1, new=FALSE)
                    plot(Nlcum[[i]]/reflength[[i]] * 100, Nl[[i]], col=i, pch=pch, xlim=xlim, ylim=ylim, xaxt="n", yaxt="n", ylab="", xlab="", main="", bty="n")  
            }
            legend("topright", legend=paste(names(N50), ": N50 = ", N50, sep=""), cex=1.2, bty="n", pch=16, pt.cex=1.2, col=seq(along=Nl)) 
            close.screen(all=TRUE)
        } 
        ## Plot cumulative contig length with ggplot2
        ## Note: ggplot2 plotting options can be looked up with theme_get()
	if(style=="ggplot2") {
                require("ggplot2")
                plotdf <- data.frame(Samples=rep(names(Nlcum), sapply(Nlcum, length)), length=unlist(Nl), perc=unlist(lapply(names(Nlcum), function(x) Nlcum[[x]]/reflength[[x]]*100)))
                counts <- table(plotdf[,1]); counts <- counts[names(N50)]
                N50rep <- paste(plotdf[,1], ": N50=", unlist(lapply(as.character(unique(plotdf[,1])), function(x) rep(N50[x], counts[names(N50[x])]))), sep="")
                plotdf[,1] <- N50rep
                ggplot(plotdf, aes(perc, length, color=Samples)) + 
                geom_point() + 
                scale_x_continuous(xlab) + 
                scale_y_continuous(ylab) + 
                labs(title = main) +
                theme(plot.title = element_text(size = sizetitle)) +
                theme(axis.title.x = element_text(size = sizex)) + 
                theme(axis.title.y = element_text(size = sizey, angle = 90)) +
                theme(legend.text = element_text(size = sizelegend)) +
                theme(legend.title = element_text(size = 0))  
        }
}

require(Biostrings)

# The first draft Italian sequence
Ital1 <- readDNAStringSet("~/Desktop/Molecular/SWD/SWD_genome/The_Italian_Jobs/drosophila_suzukii_genome_v1.fasta", format="fasta", use.names=FALSE)
head(Ital1)
str(Ital1)
(Ital1[[2,1]])
It_t <- function(rm, genome){
	if 
}



#The published Italian sequences
Ital2 <- readDNAStringSet("~/Desktop/Molecular/SWD_genome/The_Italian_Jobs/Ital_SWD_pub.fasta", format='fasta', use.names=FALSE)
head(Ital2)
str(Ital2)



#Our UCD/BGI genome
UCDv3 <- readDNAStringSet("~/Desktop/Molecular/SWD_genome/D.suzuki.contigv3/D.suzukii.scafSeq.v3", format="fasta", use.names=F)
head(UCDv3)



N <- list(assembly1=width(Ital1))
## (C) Define combined length value for each assembly set
reflength1 <- sapply(N, sum) # Example uses the combined length of the contigs for each assembly.

## (D) Compute N50 values and plot cumulative length distribution with contigStats function
pdf(file="ItalianN50.pdf", bg='white')
contigStats(N=N, reflength=reflength1, style="ggplot2") # Plots cumulative contig lengths with ggplots2 library. Alternatively, the same plot can be generated with R's base graphics by assigning "base" to the style argument. The title and axis labels can be changed with the arguments main, xlab and ylab, respectively.
dev.off()
stats <- contigStats(N=N, reflength=reflength1, style="data") # Returns a list with cumulative length vectors for each assembly and a contig summary matrix.
stats[["Contig_Stats"]] #
stats$Contig_Stats[1,2]

N2 <- list(assembly2=width(Ital2))
reflength2 <- sapply(N2, sum)
stats2 <- contigStats(N=N2, reflength=reflength2, style='data')
stats2[["Contig_Stats"]]

pdf(file="Ital_SWD_pub.pdf", bg="white")
contigStats(N=N2, reflength=reflength2, style='ggplot2')
dev.off()






#using UCD scaffold (scafSeq.v3)
N3 <- list(assembly3=width(UCDv3))
reflength3 <- sapply(N3, sum)
pdf("UCDv3.pdf", bg='white')
contigStats(N=N3, reflength=reflength3, style='ggplot2')
dev.off()

stats3 <- contigStats(N=N3, reflength=reflength3, style='data')
stats3[["Contig_Stats"]]

#for 100bp and above
UCD100 <- readDNAStringSet("~/Desktop/Molecular/SWD_genome/D.suzuki.contigv3/D.suz.scaf.100.fa", format="fasta", use.names=F)
head(UCD100)
N4 <- list(assembly4 = width(UCD100))
reflength4 <- sapply(N4, sum)
stats4 <- contigStats(N=N4, reflength=reflength4, style='data')
stats4[["Contig_Stats"]]


#for 200bp and above
UCD200 <- readDNAStringSet("~/Desktop/Molecular/SWD_genome/D.suzuki.contigv3/D.suz.scaf.200.fa", format="fasta", use.names=F)
head(UCD200)
N5 <- list(assembly5 = width(UCD200))
reflength5 <- sapply(N5, sum)
stats5 <- contigStats(N=N5, reflength=reflength5, style='data')
stats5[["Contig_Stats"]]


#500 bp cutoff
CH3 <- readDNAStringSet("~/Desktop/Molecular/SWD_genome/D.suzuki.contigv3/D.suz.scaf.500.fa", format="fasta", use.names=F)
N6 <- list(assembly6 = width(CH3))
reflength6 <- sapply(N6, sum)
stats6 <- contigStats(N=N6, reflength=reflength6, style='data')
contigStats(N=N6, reflength=reflength6, style="ggplot2")
stats6[["Contig_Stats"]]

CH4 <- readDNAStringSet("~/Desktop/Molecular/SWD_genome/D.suzuki.contigv3/D.suz.scaf.1k.fa", format="fasta", use.names=F)
N7 <- list(assembly7 = width(CH4))
reflength7 <- sapply(N7, sum)
stats7 <- contigStats(N=N7, reflength=reflength7, style='data')
contigStats(N=N7, reflength=reflength6, style="ggplot2")
stats7[["Contig_Stats"]]
hist(width(CH4))


SWD_v6 <- readDNAStringSet("~/Desktop/Molecular/SWD/SWD_genome/D.suz_v6/D.suz_v6.fa", format='fasta', use.names=F)
N8 <- list(assembly8 = width(SWD_v6))
reflength8 <- sapply(N8, sum)
stats8 <- contigStats(N=N8, reflength=reflength8, style='data')
stats8[["Contig_Stats"]]
pdf(file="SWD_v6.pdf", bg='white')
contigStats(N=N8, reflength=reflength8, style='ggplot2')
dev.off()

expectedDinucleotideFrequency <- function(x){
	bf <- alphabetFrequency(x, baseOnly=T)[DNA_BASES]
	(as.matrix(bf) %*% t(bf) - diag(bf)) / length(x)
	}
expectedDinucleotideFrequency(UCDv3)	





#SWD popgen
CH1 <- read.table("~/Desktop/Programs/popoolation_1.2.2/CH1.pi")
head(CH1)
pdf(file='Taj_pi.pdf', bg='white')
plot(CH1$V4, pch=19, ylab="Tajima's pi", las=1, xlab='Position')
abline(h=mean(CH1$V4), col='red', lwd=2, lty=2)
text(120000, 0.95, expression(paste(bar(x),"= 0.0004")), cex=1.3)
dev.off()

hist(CH1$V4, ylim=c(0, 200000), breaks=20, col='grey', xlab="Tajima's pi", main="")
mean(CH1$V4)



CH1a <- read.table("~/Desktop/Programs/popoolation_1.2.2/CH1_2.pi")
mean(CH1a$V4)



CH2 <- read.table("~/Desktop/Programs/popoolation_1.2.2/CHI_sub.TajD")
head(CH2)

plot(CH2$V4, pch=19, ylab="Tajima's D", las=1, xlab='Position', ylim=c(-1,1))

#SWD mito China
CH_m <- read.table("~/Desktop/Molecular/SWD_genome/CHI_SWD_mito.pi")
head(CH_m)
length(CH_m$V4)

pdf(file="CH_mt_pi.pdf", bg='white')
plot(1:150, seq(0,1, by=0.006666667), type="n", xaxt="n", ylab="Tajima's pi", xlab="mtDNA position", las=1, cex.lab=1.2, main="China")
axis(1, at=c(1, 25, 50, 75, 100, 125, 150), lab=expression(0, "25K", "50K", "75K", "100K", "125K", "150K"), cex.axis=1.2, cex.lab=1.2)
points(CH_m$V4, pch=19)
abline(h=mean(CH_m$V4), col='red', lwd=3, lty=2)
dev.off()





#SWD Wolbachia
SWD_wol<- readDNAStringSet("~/Desktop/Molecular/SWD_Wolbachia/Dsuz_wol.fa", format="fasta", use.names=F)
head(SWD_wol)
N7 <- list(assembly7 = width(SWD_wol))
reflength7 <- sapply(N7, sum)
stats7 <- contigStats(N=N7, reflength=reflength7, style='data')
contigStats(N=N7, reflength=reflength7, style="ggplot2")
stats7[["Contig_Stats"]]







#China wolbachia, Tajima's pi
CH_w1 <- read.table("~/Desktop/Molecular/SWD_Wolbachia/CHI_wol1.pi")
head(CH_w1)
plot(CH_w1$V4)#wow! no way, I don't believe this
mean(CH_w1$V4)
plot(CH_w1$V5)


#for min 2 snps per window
CH_w2 <- read.table("~/Desktop/Molecular/SWD_Wolbachia/CHI_wol2.pi")
head(CH_w2)
plot(CH_w2$V4)
mean(CH_w2$V4)#must be shitty mapping or too low coverage

#require 20x covereage and min 10 snps per window (1k window)
CH_w10 <- read.table("~/Desktop/Molecular/SWD_Wolbachia/CHI_wol10.pi")
head(CH_w10)
plot(CH_w10$V4)#still super high. 

#require 20x coverage, 10 snps per window (100bp window)
CH_w10b <- read.table("~/Desktop/Molecular/SWD_Wolbachia/CHI_wol10b.pi")
head(CH_w10b)
plot(CH_w10b$V4)
mean(CH_w10b$V4)


#CHI_6 Tajima's pi
CHI_6 <- read.table("~/Desktop/CHI_6.pi")#pool size 50
head(CHI_6)
par(mfrow=c(2,2))
plot(CHI_6$V4, ylab="Tajima's pi", las=1, pch=19, main="10k win mincnt1 mincov4 maxcov1k", ylim=c(0, 0.06))
mean(CHI_6$V4)
length(CHI_6$V4)

CHI_6_min2 <- read.table("~/Desktop/CHI_6_min2.pi")#pool size 50
head(CHI_6_min2)
length(CHI_6_min2$V4)
plot(CHI_6_min2$V4, ylab="Tajima's pi", las=1, pch=19, main='10win mincnt1 mincov2 maxcov1K', ylim=c(0, 0.06))
mean(CHI_6_min2$V4)







#Not sure about this one. Scale is funny
CHI_6_1k <- read.table("~/Desktop/CHI_6_1k.pi")#pool size 50
head(CHI_6_1k)
length(CHI_6_1k$V4)
plot(CHI_6_1k$V4, ylab="Tajima's ", las=1, pch=19, main="1kwin mincnt1 mincov2 maxcov10k",ylim=c(0, 0.6))#window size effects scale of pi but not mean?
mean(CHI_6_1k$V4)

#Increasing poolsize to 100, should be number of c'somes
CHI_6_100 <- read.table("~/Desktop/CHI_6_100.pi")
head(CHI_6_100)
length(CHI_6_100$V4)
plot(CHI_6_100$V4, ylab="Tajima's Pi", las=1, pch=19, ylim=c(0, 0.06))
mean(CHI_6_100$V4)

#Tajima's D for China samples
CHI_6_D <- read.table("~/Desktop/CHI_6_1k.D")
head(CHI_6_D)
plot(CHI_6_D$V4, ylim=c(-0.1,0.1), pch=19, ylab="Tajima's D")





#Tajima's pi for WAT and CHI with bowtie 2

CHI_1.pi <- read.table("~/Desktop/Molecular/SWD/SWD_genome/CHI_1.pi")
head(CHI_1.pi)
length(CHI_1.pi$V5)

range(CHI_1.pi$V5)
plot(1:length(CHI_1.pi$V5), CHI_1.pi$V5, ylab="Tajima's pi", las=1, type="line")
abline(h=mean(CHI_1.pi$V5), col="red") #0.669

CHI_1.pi$V5 <- as.numeric(as.character(CHI_1.pi$V5))
 
pdf(file="Chinapi.pdf", bg='white')
plot(CHI_1.pi$V5, ylab="Tajima's pi", las=1, pch=19, cex=0.2, main="China", ylim=c(0, 0.07))
mean(CHI_1.pi$V5, na.rm=T)
abline(h=median(CHI_1.pi$V5, na.rm=T), col='red', lwd=2)
dev.off()

plot(CHI_1.pi$V5, ylab="Tajima's pi", las=1, type='l', main="China", ylim=c(0, 0.07))



CHI_1A.pi <- read.table("~/Desktop/Molecular/SWD/SWD_genome/CHI_1A.pi")
head(CHI_1A.pi)
length(CHI_1A.pi$V4)
plot(CHI_1A.pi$V4, ylab="Tajima's pi", las=1, type="line")
abline(h=mean(CHI_1A.pi$V4), col="red") #0.607
plot(CHI_1A.pi$V4, ylab="Tajima's pi", las=1, pch=19, cex=0.2)

CHI_1b.pi <- read.table("~/Desktop/Molecular/SWD/SWD_genome/CHI_1b.pi")
head(CHI_1b.pi)
length(CHI_1b.pi$V4)
plot(CHI_1b.pi$V4, ylab="Tajima's pi", las=1, type="line")
mean(CHI_1b.pi$V4)
abline(h=mean(CHI_1b.pi$V4), col="red") #0.536
plot(CHI_1b.pi$V4, ylab="Tajima's pi", las=1, pch=19, cex=0.2)

CHI_1c.pi <- read.table("~/Desktop/Molecular/SWD/SWD_genome/CHI_1c.pi")
head(CHI_1c.pi)
length(CHI_1c.pi$V4)
plot(CHI_1c.pi$V4, ylab="Tajima's pi", las=1, type="line")
mean(CHI_1c.pi$V4)
abline(h=mean(CHI_1c.pi$V4), col="red")
plot(CHI_1c.pi$V4, ylab="Tajima's pi", las=1, pch=19, cex=0.2)

plot(CHI_1c.pi$V3/CHI_1c.pi$V2, cex=0.2, pch=19, ylab='', las=1)
mean(CHI_1c.pi$V3/CHI_1c.pi$V2)
abline(h=mean(CHI_1c.pi$V3/CHI_1c.pi$V2), lwd=2, col='red')
#WAT

WAT_1.pi <- read.table("~/Desktop/Molecular/SWD/SWD_genome/WAT1_1.pi")
head(WAT_1.pi)

WAT_1.pi$V5 <- as.numeric(as.character(WAT_1.pi$V5))

pdf(file="Watpi.pdf", bg='white')
plot(WAT_1.pi$V5, las=1, ylab="Tajima's pi", main="Watsonvile", pch=19, cex=0.2, ylim=c(0, 0.07))
mean(WAT_1.pi$V5, na.rm=T)#0.013
abline(h=median(WAT_1.pi$V5, na.rm=T), col="red", lwd=2)
dev.off()

plot(WAT_1.pi$V5, las=1, ylab="Tajima's pi", main="Watsonvile", type='l', ylim=c(0, 0.07))


length(WAT_1.pi$V5)
par(mfrow=c(1,2))
plot(WAT_1.pi$V5, ylab="Tajima's pi", las=1, type="line")
abline(h=mean(WAT_1.pi$V5, na.rm=R), col="red") #0.567

#need to account for window size 
WAT_1a.pi <- read.table("~/Desktop/Molecular/SWD/SWD_genome/WAT1_1a.pi")
head(WAT_1a.pi)
length(WAT_1a.pi$V4)
par(mfrow=c(1,2))
plot(WAT_1a.pi$V4, ylab="Tajima's pi", las=1, type="line")
abline(h=mean(WAT_1a.pi$V4), col="red") #0.47
plot(WAT_1a.pi$V4, ylab="Tajima's pi", las=1, pch=19, cex=0.2)

WAT_1b.pi <- read.table("~/Desktop/Molecular/SWD/SWD_genome/WAT1_1b.pi")
head(WAT_1b.pi)
length(WAT_1b.pi$V4)
par(mfrow=c(1,2))
plot(WAT_1b.pi$V4, ylab="Tajima's pi", las=1, type="line")
mean(WAT_1b.pi$V4)
abline(h=mean(WAT_1b.pi$V4), col="red") #0.47
plot(WAT_1b.pi$V4, ylab="Tajima's pi", las=1, pch=19, cex=0.2)
abline(h=mean(WAT_1b.pi$V4), col="red")

WAT_1c.pi <- read.table("~/Desktop/Molecular/SWD/SWD_genome/WAT1_1c.pi")
head(WAT_1c.pi)
length(WAT_1c.pi$V4)
par(mfrow=c(1,2))
plot(WAT_1c.pi$V4, ylab="Tajima's pi", las=1, type="line")
mean(WAT_1c.pi$V4)
abline(h=mean(WAT_1c.pi$V4), col="red")
plot(WAT_1c.pi$V4, ylab="Tajima's pi", las=1, pch=19, cex=0.2)
abline(h=mean(WAT_1c.pi$V4), col="red") 



#Tajima's D with initial settings
#CHI_1
CHI_1.D <- read.table("~/Desktop/Molecular/SWD/SWD_genome/CHI_1.D")
head(CHI_1.D)
length(CHI_1.D$V5)

CHI_1.D$V5 <- as.numeric(as.character(CHI_1.D$V5))
median(CHI_1.D$V5, na.rm=T)


pdf(file='ChiD.pdf', bg='white')
plot(CHI_1.D$V5, las=1, ylab="Tajima's D", pch=19, cex=0.2, ylim=c(-3, 1), main="China")
mean(CHI_1.D$V5, na.rm=T)-0.575
abline(h=median(CHI_1.D$V5, na.rm=T), col='red', lwd=2)
dev.off()

#WAT_1
WAT_1.D <- read.table("~/Desktop/Molecular/SWD/SWD_genome/WAT1_1.D")
head(WAT_1.D)
length(WAT_1.D$V5)

WAT_1.D$V5 <- as.numeric(as.character(WAT_1.D$V5))
length(WAT_1.D$V5)

median(WAT_1.D$V5, na.rm=T)

pdf(file='WatD.pdf', bg='white')
plot(WAT_1.D$V5, ylab="Tajima's D", las=1, pch=19, cex=0.2, ylim=c(-3, 1), main="Watsonville")
mean(WAT_1.D$V5, na.rm=T)
abline(h=median(WAT_1.D$V5, na.rm=T), col="red", lwd=2) #-0.400
dev.off()



#Wat-Chi Fst from v3 genome
Fst1 <- read.delim("~/Desktop/WAT_CHI_fst.txt", header=F)
str(Fst1)
head(Fst1)

Fst1 <- as.numeric(as.character(Fst1$V1))
mean(Fst1, na.rm=T)
median(Fst1, na.rm=T)

pdf(file="Wat_CHI_fst.pdf", bg="white")
plot(Fst1, ylab="Fst", las=1, cex=0.2)
abline(h=median(Fst1, na.rm=T), lwd=2, col="red")
dev.off()


#pairwaise allele freqs
pwc1 <- read.delim("~/Desktop/WAT.CHI_1_pwc", header=T)
str(pwc1)

pwc1 <- as.numeric(as.character(pwc1$diff.1.2))
summary(pwc1)

plot(pwc1, las=1, cex=0.2, ylab='Allele freqs')




setwd("~/Desktop/Molecular/SWD/SWD_genome/SWD_ortho/5203")
e1 <- read.table("codon.coa", header=T, sep=" ")
pdf(file='5203_eigen.pdf', bg='white')
plot(e1[,2], e1[,3], pch=19, ylab='Axis 2', xlab='Axis 1', las=1)
abline(h=0, v=0, lwd=1)
dev.off()
