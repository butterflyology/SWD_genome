#SWD ortholog analysis
library(plyr)
library(seqinr)
library(exactRankTests)
library(ade4)
library(vcd)
setwd("~/Desktop/Projects/SWD_genome/SWD_genome_data/SWD_genome")
#the CAI for this is not appropriate, based on a bacterium
load("Dros_dnds.R")

#####Orthologs
Ortho <- read.delim("SWD_ortho_out.txt", sep="\t")
Ortho$Nc <- as.numeric(Ortho$Nc)#something is wrong with Nc, many values not divided by 64 (apparently)

W <- Ortho$Nc
Nc=NA
for(t in 1:length(W)){
	if(W[t] > 64) Nc[t] = (W[t]/64)
	else Nc[t] = Nc[t]
} 

Ortho$Nc <- Nc
str(Ortho)
head(Ortho)
#plot(Ortho)
par(mfrow=c(1,2))
hist(Ortho$Nc, col='grey', xlim=c(0, 64), las=1)

pdf(file="L_aa.pdf", bg="white")
hist(Ortho$L_aa, col="grey", xlab="AA length", main="", xlim=c(0,10000), las=1)
abline(v=median(Ortho$L_aa), col="red", lwd=2, lty=2)
legend("topright", legend=c("Median length"), lwd=2, lty=2, bty="n", col="red")
dev.off()

Dsuz <- subset(Ortho, Taxon=="Dsuz")
Dsuz <- droplevels(Dsuz)
Dbia <- subset(Ortho, Taxon == "Dbia")
Dbia <- droplevels(Dbia)
Dtak <- subset(Ortho, Taxon == "Dtak")
Dtak <- droplevels(Dtak)
Dmel <- subset(Ortho, Taxon == "Dmel")
Dmel <- droplevels(Dmel)


#######FPKM
fpkm <- read.csv("Ortho_out.csv")
str(fpkm)
head(fpkm)
sum(Dsuz$Gene %in% fpkm$Gene)

#merge Dsuz and fpkm
Dsuz_fpkm <- merge(Dsuz, fpkm, by="Gene")
str(Dsuz_fpkm)
head(Dsuz_fpkm)
length(Dsuz_fpkm$Gene)
summary(Dsuz_fpkm$G3s)
summary(Dsuz_fpkm$Male)
summary(Dsuz_fpkm$Female)


######## X-linked genes
xgenes <- read.delim("~/Desktop/Programs/ncbi-blast-2.2.27+/databases/Dmel/SWD_u.txt", sep="\t")
str(xgenes)
head(xgenes)
length(unique(xgenes$Gene))
Uniq_X <- as.matrix(unique(xgenes$Gene))
colnames(Uniq_X) <- "Gene"
#intersect(Ortho$Gene, Uniq_X)
sum(Ortho$Gene %in% xgenes$Gene)


Dros_X <- merge(Ortho, Uniq_X, by="Gene")
Dros_X <- droplevels(Dros_X)
str(Dros_X)
head(Dros_X)
tail(Dros_X)

#write.table(Dros_X, file="Dros_x.txt", sep="\t", eol="\r", quote=FALSE)

#subset X-linked genes by taxon
Dsuz_X <- droplevels(subset(Dros_X, Taxon=="Dsuz"))
str(Dsuz_X)
length(Dsuz_X$Gene) #802 X-linked genes
Dbia_X <- droplevels(subset(Dros_X, Taxon=="Dbia"))
str(Dbia_X) #802 X-linked genes
Dtak_X <- droplevels(subset(Dros_X, Taxon=="Dtak"))
str(Dtak_X) #801 X-linked genes
Dmel_X <- droplevels(subset(Dros_X, Taxon=="Dmel"))
str(Dmel_X) #937 X-linked genes


######## Autosomal genes
Auto <- as.matrix(setdiff(Ortho$Gene, xgenes$Gene))
colnames(Auto) <- "Gene"
Dros_A <- merge(Ortho, Auto, by="Gene")
str(Dros_A)
head(Dros_A)
#write.table(Dros_A, file="Dros_A.txt", sep="\t", eol="\r", quote=FALSE)


#Subset by taxon
Dsuz_A <- droplevels(subset(Dros_A, Taxon=="Dsuz"))
str(Dsuz_A) #4117 Autosomal genes
Dbia_A <- droplevels(subset(Dros_A, Taxon=="Dbia"))
str(Dbia_A) #4117 Autosomal genes
Dtak_A <- droplevels(subset(Dros_A, Taxon=="Dtak"))
str(Dtak_A) #4118 Autosomal genes
Dmel_A <- droplevels(subset(Dros_A, Taxon=="Dmel"))
str(Dmel_A) #3982 Autosomal genes

#Don't believe the CAI
#CBI/Fop set to Drosophila melanogaster
hist(Ortho$CBI, col="grey", xlim=c(-0.5, 1.0), main="", las=1, xlab=expression(paste(italic("Drosophila"), " CBI")))
hist(Ortho$GC3, col='grey', xlim=c(0,1), ylim=c(0, 4000), main="", las=1, xlab=expression(paste(italic("Drosophila")," GC content")))

boxplot(Dsuz$CBI, Dbia$CBI, Dtak$CBI, Dmel$CBI, names=c("D. suz", "D. bia", "D. tak", "D. mel"), col="grey", pch=19, cex=0.5, las=1, ylim=c(-0.5, 1.0), main=expression(paste(italic("Drosophila"), " CBI")))

pdf(file="A_X_CBI.pdf", bg="white")
par(mfrow=c(1,2))
boxplot(Dsuz_A$CBI, Dbia_A$CBI, Dtak_A$CBI, Dmel_A$CBI, names=c("D. suz", "D. bia", "D. tak", "D. mel"), col="grey", pch=19, cex=0.5, las=1, ylim=c(-0.5, 1.0), main=expression(paste(italic("Drosophila"), " autosomal CBI")))

boxplot(Dsuz_X$CBI, Dbia_X$CBI, Dtak_X$CBI, Dmel_X$CBI, names=c("D. suz", "D. bia", "D. tak", "D. mel"), col="grey", pch=19, cex=0.5, las=1, ylim=c(-0.5, 1.0), main=expression(paste(italic("Drosophila"), " X CBI")))
dev.off()


wilcox.exact(Dsuz_A$CBI, Dbia_A$CBI, conf.int=TRUE, conf.level=0.95, na.rm=TRUE, exact=FALSE)

wilcox.exact(Dsuz_X$CBI, Dbia_X$CBI, conf.int=TRUE, conf.level=0.95, na.rm=TRUE, exact=FALSE)




boxplot(Dsuz$GC3s, Dbia$GC3s, Dtak$GC3s, Dmel$GC3s, names=c("D. suz", "D. bia", "D. tak", "D. mel"), ylim=c(0,1), col="grey", pch=19, cex=0.5, las=1, main=expression(paste(italic("Drosophila"), " GC3"[s])))
abline(h=median(Dmel$GC3), lwd=2, lty=2, col="red")
t.test(Dsuz$GC3s, Dbia$GC3s)



boxplot(Dsuz$GC, Dbia$GC, Dtak$GC, Dmel$GC, names=c("D. suz", "D. bia", "D. tak", "D. mel"), col="grey", pch=19, cex=0.5, ylim=c(0,1), las=1, main=expression(paste(italic("Drosophila"), " GC")))
abline(h=median(Dbia$GC), lwd=2, lty=2, col="red")

boxplot(Dsuz$Fop, Dbia$Fop, Dtak$Fop, Dmel$Fop, names=c("Dsuz", "Dbia", "Dtak", "Dmel"), ylim=c(0,1), col="grey", pch=19, cex=0.5)
abline(h=median(Dbia$Fop), lwd=2, lty=2, col="red")

boxplot(Dsuz$Nc, Dbia$Nc, Dtak$Nc, Dmel$Nc, names=c("D. suz", "D. bia", "D. tak", "D. mel"), col="grey", pch=19, cex=0.8, ylim=c(0, 64), main=expression(paste(italic("Drosophila"), " eNc")))


plot(Dsuz$GC, Dsuz$GC3, xlab="GC content", ylab=expression(paste("GC"[3], " content")), pch=19, cex=0.5, ylim=c(0.15, 1), xlim=c(0.2,0.8))
(lm1 <- lm(Dsuz$GC3~Dsuz$GC))
abline(lm1, lwd=2, lty=2, col='red')







plot(log10(fpkm$Female), log10(fpkm$Male), xlim=c(-4, 4), ylim=c(-4,4), xlab=expression(paste("log"["10" ,] ," female " ,italic("D. suzukii"))), ylab=expression(paste("log"["10" ,] ," male " ,italic("D. suzukii"))), type="n")
points(log10(fpkm$Female), log10(fpkm$Male), pch=19, cex=0.8, col=c("red", "blue"))
abline()


(lm(log10(fpkm$Female) ~ log10(fpkm$Male), na.action=))




plot(log10(Dsuz_fpkm$Male), Dsuz_fpkm$G3s, pch=19, cex=0.8, xlab=expression(paste("log"["10" ,] ," male " ,italic("D. suzukii"))), ylab=expression(paste(italic("D. suzukii "), " GC"["3"] )), las=1)
#text(-4, 0.6, "\\VE", vfont=c("sans serif symbol", "plain"), cex=2)
summary((lm(log10(Dsuz_fpkm$G3s) ~ Dsuz_fpkm$Male)))
summary((lm(Dsuz_fpkm$G3s ~ Dsuz_fpkm$Female)))

plot(log10(Dsuz_fpkm$Female), Dsuz_fpkm$G3s, pch=19, cex=0.8, xlab=expression(paste("log"["10" ,] ," female " ,italic("D. suzukii"))), ylab=expression(paste(italic("D. suzukii "), " GC"["3"] )), las=1)

par(mfrow=c(1,2))
hist(log10(Dsuz_fpkm$Female), col='grey', xlim=c(-6, 6), main="", las=1, xlab=expression(paste("log"["10" ,] ," female " ,italic("D. suzukii"))))
hist(log10(Dsuz_fpkm$Male), col='grey', xlim=c(-6, 6), main="", las=1, xlab=expression(paste("log"["10" ,] ," male " ,italic("D. suzukii"))))

#pdf(file="GC3_fpkm.pdf", bg="white")
par(mfrow=c(1,2))
plot(log(Dsuz_fpkm$Male), Dsuz_fpkm$GC3, ylim=c(0,1), pch=19, las=1, cex=0.8, ylab=expression("GC"[3]), xlab=expression(paste("log"["10" ,], " male FPKM")))
text(-11, 1.0, "A)")
plot(log(Dsuz_fpkm$Female), Dsuz_fpkm$GC3, ylim=c(0,1), pch=19, las=1, cex=0.8, ylab="", xlab=expression(paste("log"["10" ,], " female FPKM")), yaxt="n")
text(-4.8, 1.0, "B)")
#dev.off()



boxplot(Dros_X$GC3s, Dros_A$GC3s, names=c("X", "Autosome"), ylim=c(0,1), col="grey", pch=19, cex=0.8, main=expression(paste(italic("Drosophila"), " GC"[3])), las=1)

par(mfrow=c(1,2))
#pdf(file="Dsuz_GC3.pdf", bg="white")
boxplot(Dsuz_X$GC3s, Dsuz_A$GC3s, names=c("X", "Autosome"), ylim=c(0,1), col='grey', pch=19, cex=0.8, main=expression(paste(italic("D. suzukii"))), ylab=expression(paste("GC"[3])), las=1)
#dev.off()

plot(Dsuz$GC3s, Dsuz$CBI, pch=19, xlab=expression(paste("GC3"[s])), ylab="CBI", las=1)

wilcox.exact(Dsuz$GC3s, Dsuz$CBI, conf.int=TRUE, conf.level=0.95, na.rm=TRUE, exact=FALSE)
(lm1 <- lm(Dsuz$GC3s~Dsuz$CBI))
summary(lm1)

boxplot(Dmel_X$GC3s, Dmel_A$GC3s, names=c("X", "Autosome"), ylim=c(0,1), col='grey', pch=19, cex=0.8, main=expression(paste(italic("D. melanogaster"), " GC"[3])))

boxplot(Dsuz_X$GC3, Dbia_X$GC3, Dtak_X$GC3, Dmel_X$GC3, cex=0.8, pch=19, col="grey")


wilcox.exact(Dsuz_X$GC3s, Dsuz_A$GC3s, conf.int=TRUE, conf.level=0.95, na.rm=TRUE, exact=FALSE)


hist(Dsuz_A$GC3s, xlim=c(0,1), col="grey", las=1, xlab=expression(paste(italic( "D. suzukii"), " GC"[3])), main="")
abline(v=median(Dsuz_A$GC3s), lty=2, lwd=2)
hist(Dsuz_X$GC3s, xlim=c(0,1), col="red", add=TRUE)
abline(v=median(Dsuz_X$GC3s), lty=2, lwd=2, col="red")
legend("topleft", legend=c("Autosome", "X"), bty="n", lty=2, lwd=2, col=c("black", "red"))



perm.test(Dsuz$Fop, Dbia$Fop, conf.int=TRUE, conf.level=0.95, na.rm=TRUE, exact=FALSE)

boxplot(Dsuz_X$Fop, Dbia_X$Fop, names=c("Dsuz X Fop", "Dbia X Fop"), ylim=c(0,1), col='grey', pch=19, cex=0.8)
perm.test(Dsuz_X$Fop, Dbia_X$Fop, conf.int=TRUE, conf.level=0.95, na.rm=TRUE, exact=FALSE)
wilcox.exact(Dsuz_X$Fop, Dbia_X$Fop, conf.int=TRUE, conf.level=0.95, na.rm=TRUE)

boxplot(Dsuz_A$Fop, Dbia_A$Fop, names=c("Dsuz A Fop", "Dbia A Fop"), ylim=c(0,1), col='grey', pch=19, cex=0.8)
perm.test(Dsuz_A$Fop, Dbia_A$Fop, conf.int=TRUE, conf.level=0.95, na.rm=TRUE, exact=FALSE)
wilcox.exact(Dsuz_A$Fop, Dbia_A$Fop, conf.int=TRUE, conf.level=0.95, na.rm=TRUE)


#####################
#calculating dN/dS
require(seqinr)
require(plyr)
setwd("~/Desktop/Molecular/SWD/SWD_genome/SWD_BEAST/fasta")
g2300 <- read.alignment("FAMILY.2300.aligned_codon", format="fasta")

k2300 <- kaks(g2300)
k2300
k1 <- k2300$ka/k2300$ks
str(k1)
k1

files <- list.files(path="~/Desktop/Molecular/SWD/SWD_genome/SWD_BEAST/fasta", pattern="*aligned_codon", full.names=T, recursive=FALSE)

mas1 <- laply(files, function(x){
	t <- read.alignment(x, format="fasta")
	out <- kaks(t)
	ka_ks <- out$ka/out$ks
	mka_ks <- as.matrix(ka_ks)
	mka_ks <- mka_ks[,c("Dsuz", "Dbia", "Dtak", "Dmel")]
	mka_ks <- mka_ks[c("Dsuz", "Dbia", "Dtak", "Dmel"),]
	return(mka_ks)
	#write.table(mka_ks, "~/Desktop/Molecular/SWD/SWD_genome/SWD_BEAST/fasta/SWD_kaks.txt", append=TRUE, row.names=TRUE, col.names=TRUE, , eol="\r")
})

which(mas1 == Inf)
which(mas1 < 0)
mas1[which(mas1 == Inf)] <- NA
mas1[which(mas1 < 0)] <- NA
aaply(mas1, c(2,3), function(x) mean(x, na.rm=TRUE))

require(rethinking)
dens(mas1[,2,1], las=1, )
#array dims are: dimension of list (4919 here), columns, row
#[,2,1] is Dtak Dbia
boxplot(mas1[,2,1], mas1[,3,1], mas1[,4,1], names=c("Dsuz-Dbia", "Dsuz-Dtak", "Dsuz-Dmel"), las=1, col="grey", pch=19, cex=0.8, ylab="dN / dS" )

#triangle plots


#just dN and DS
m_a <- laply(files, function(x){
	t <- read.alignment(x, format="fasta")
	out <- kaks(t)
	out_ka <- out$ka
	out_ka <- as.matrix(out_ka)
	out_ka <- out_ka[,c("Dsuz", "Dbia", "Dtak", "Dmel")]
	out_ka <- out_ka[c("Dsuz", "Dbia", "Dtak", "Dmel"),]
	return(out_ka)
})

m_a[which(m_a == Inf)] <- NA
m_a[which(m_a < 0)] <- NA
m_a[which(m_a > 9)] <- NA
aaply(m_a, c(2,3), function(x) mean(x, na.rm=TRUE))

#autosomal gene ka
mcols <- data.frame(cbind(m_a[,1,2], m_a[,2,3], m_a[,3,1]))
mcols <- mcols[complete.cases(mcols),]
colnames(mcols) <- c("suz - bia", "bia - tak", "tak - suz")

#remove row when all columns contain a 0
mcols <- mcols[apply(mcols[c(1:3)], 1, function(z) any(z!=0)),]

pdf(file="all_genes_dn.pdf", bg="white")
triangle.plot(mcols, cpoint=0.8, show.position=FALSE, scale=FALSE, min3=c(0,0,0), max3=c(1,1,1))
dev.off()



m_s <- laply(files, function(x){
	t <- read.alignment(x, format="fasta")
	out <- kaks(t)
	out_ks <- out$ks
	out_ks <- as.matrix(out_ks)
	out_ks <- out_ks[,c("Dsuz", "Dbia", "Dtak", "Dmel")]
	out_ks <- out_ks[c("Dsuz", "Dbia", "Dtak", "Dmel"),]
	return(out_ks)
})

m_s[which(m_s == Inf)] <- NA
m_s[which(m_s < 0)] <- NA
m_s[which(m_s > 9)] <- NA

aaply(m_s, c(2,3), function(x) mean(x, na.rm=TRUE))

boxplot(m_a[,2,1], m_a[,3,1], m_a[,4,1],las=1, pch=19, cex=0.8, col="grey", names=c("Dsuz - Dbia", "Dsuz - Dtak", "Dsuz - Dmel"), ylab= "A-chromosome dN")

boxplot(m_s[,2,1], m_s[,3,1], m_s[,4,1],las=1, pch=19, cex=0.8, col="grey", names=c("Dsuz - Dbia", "Dsuz - Dtak", "Dsuz - Dmel"), ylab= "A-chromosome dS")

#Take tope 1% of ka genes (~50)
m_a[which(m_a > 0.67)]



#break down by autosomal and X

#X-linked genes

file1 <- list.files(path="~/Desktop/Molecular/SWD/SWD_genome/Xout", pattern="*aligned_codon", full.names=T, recursive=FALSE)

xdnds <- laply(file1, function(x){
	t <- read.alignment(x, format="fasta")
	out <- kaks(t)
	ka_ks <- out$ka/out$ks
	xka_ks <- as.matrix(ka_ks)
	xka_ks <- xka_ks[,c("Dsuz", "Dbia", "Dtak", "Dmel")]
	xka_ks <- xka_ks[c("Dsuz", "Dbia", "Dtak", "Dmel"),]
	return(xka_ks)
})

head[which(xdnds == Inf)] <- NA
xdnds[which(xdnds < 0)] <- NA
aaply(xdnds, c(2,3), function(x) median(x, na.rm=TRUE))


xdn <- laply(file1, function(x){
	t <- read.alignment(x, format="fasta")
	out <- kaks(t)
	ka <- out$ka
	xka <- as.matrix(ka)
	xka <- xka[,c("Dsuz", "Dbia", "Dtak", "Dmel")]
	xka <- xka[c("Dsuz", "Dbia", "Dtak", "Dmel"),]
	return(xka)
})
xdn[which(xdn == Inf)] <- NA
xdn[which(xdn < 0)] <- NA
aaply(xdn, c(2,3), function(x) median(x, na.rm=TRUE))


boxplot(xdn[,2,1], xdn[,3,1], xdn[,3,2], xdn[,4,1],las=1, pch=19, cex=0.8, col="grey", names=c("Dsuz - Dbia", "Dsuz - Dtak", "Dbia-Dtak", "Dsuz - Dmel"), ylab= "X- dN")

xdn[which(xdn > 0.8)]#874




xds <- laply(file1, function(x){
	t <- read.alignment(x, format="fasta")
	out <- kaks(t)
	ks <- out$ks
	xks <- as.matrix(ks)
	xks <- xks[,c("Dsuz", "Dbia", "Dtak", "Dmel")]
	xks <- xks[c("Dsuz", "Dbia", "Dtak", "Dmel"),]
	return(xks)
})
xds[which(xds == Inf)] <- NA
xds[which(xds < 0)] <- NA
xds[which(xds > 9)] <- NA
aaply(xds, c(2,3), function(x) median(x, na.rm=TRUE))

boxplot(xds[,2,1], xds[,3,1], xds[,3,2], xds[,4,1],las=1, pch=19, cex=0.8, col="grey", names=c("Dsuz - Dbia", "Dsuz - Dtak", "Dbia-Dtak", "Dsuz - Dmel"), ylab= "X- dS", ylim=c(0, 2))

xcols <- data.frame(cbind(xdn[,1,2], xdn[,2,3], xdn[,3,1]))
xcols <- xcols[complete.cases(xcols),]
colnames(xcols) <- c("suz - bia", "bia - tak", "tak - suz")
#remove row when all columns contain a 0
xcols <- xcols[apply(xcols[c(1:3)], 1, function(z) any(z!=0)),]

pdf(file="x_genes_dn.pdf", bg="white")
triangle.plot(xcols, cpoint=0.8, show.position=FALSE, label=namies, scale=FALSE, min3=c(0,0,0), max3=c(1,1,1))
dev.off()





########pull top 1% genes for dNs for whole genome and look at what they are########

boxplot(xdnds[,2,1], xdnds[,3,1], xdnds[,3,2], xdnds[,4,1],las=1, ylim=c(0,1), pch=19, cex=0.8, col="grey", names=c("Dsuz - Dbia", "Dsuz - Dtak", "Dbia-Dtak", "Dsuz - Dmel"), ylab= "X-chromosome dN/dS")


#Autosomal data
file2 <- list.files(path="~/Desktop/Molecular/SWD/SWD_genome/Aout", pattern="*aligned_codon", full.names=T, recursive=FALSE)

adnds <- laply(file2, function(x){
	t <- read.alignment(x, format="fasta")
	out <- kaks(t)
	ka_ks <- out$ka/out$ks
	aka_ks <- as.matrix(ka_ks)
	aka_ks <- aka_ks[,c("Dsuz", "Dbia", "Dtak", "Dmel")]
	aka_ks <- aka_ks[c("Dsuz", "Dbia", "Dtak", "Dmel"),]
	return(aka_ks)
})

adnds[which(adnds == Inf)] <- NA
adnds[which(adnds < 0)] <- NA
aaply(adnds, c(2,3), function(x) median(x, na.rm=TRUE))

boxplot(adnds[,2,1], adnds[,3,1], adnds[,3,2], adnds[,4,1],las=1, pch=19, cex=0.8, col="grey", names=c("Dsuz - Dbia", "Dsuz - Dtak", "Dbia-Dtak", "Dsuz - Dmel"), ylab= "Autosomal dN/dS")



adn <- laply(file2, function(x){
	t <- read.alignment(x, format="fasta")
	out <- kaks(t)
	ka <- out$ka
	aka <- as.matrix(ka)
	aka <- aka[,c("Dsuz", "Dbia", "Dtak", "Dmel")]
	aka <- aka[c("Dsuz", "Dbia", "Dtak", "Dmel"),]
	return(aka)
})
adn[which(adn == Inf)] <- NA
adn[which(adn < 0)] <- NA
aaply(adn, c(2,3), function(x) median(x, na.rm=TRUE))

boxplot(adn[,2,1], adn[,3,1], adn[,3,2], adn[,4,1],las=1, pch=19, cex=0.8, col="grey", names=c("Dsuz - Dbia", "Dsuz - Dtak", "Dbia-Dtak", "Dsuz - Dmel"), ylab= "Autosomal dN")

acols <- data.frame(cbind(adn[,1,2], adn[,2,3], adn[,3,1]))
acols <- acols[complete.cases(acols),]
colnames(acols) <- c("suz - bia", "bia - tak", "tak - suz")

#remove row when all columns contain a 0
acols <- acols[apply(acols[c(1:3)], 1, function(z) any(z!=0)),]

pdf(file="auto_gene_dn.pdf", bg="white")
triangle.plot(acols, cpoint=0.8, show.position=FALSE, scale=FALSE, min3=c(0,0,0), max3=c(1,1,1))
dev.off()

pdf(file="auto_x_dn.pdf", bg="white")
par(mfrow=c(2,1))

#autosomal
triangle.plot(acols, cpoint=0.8, show.position=FALSE, scale=FALSE, min3=c(0,0,0), max3=c(1,1,1))

#x-linked
triangle.plot(xcols, cpoint=0.8, show.position=FALSE, label=namies, scale=FALSE, min3=c(0,0,0), max3=c(1,1,1))
dev.off()





ads <- laply(file2, function(x){
	t <- read.alignment(x, format="fasta")
	out <- kaks(t)
	ks <- out$ks
	aks <- as.matrix(ks)
	aks <- aks[,c("Dsuz", "Dbia", "Dtak", "Dmel")]
	aks <- aks[c("Dsuz", "Dbia", "Dtak", "Dmel"),]
	return(aks)
})
ads[which(ads == Inf)] <- NA
ads[which(ads < 0)] <- NA
ads[which(ads > 8)] <- NA
aaply(ads, c(2,3), function(x) median(x, na.rm=TRUE))

boxplot(ads[,2,1], ads[,3,1], ads[,3,2], ads[,4,1],las=1, pch=19, cex=0.8, col="grey", names=c("Dsuz - Dbia", "Dsuz - Dtak", "Dbia-Dtak", "Dsuz - Dmel"), ylab= "Autosomal dS")



save(file="Dros_dnds.R", list=c("adnds", "xdnds", "m_s", "m_a", "mas1", "Ortho", "Dsuz", "Dmel", "Dbia", "Dtak", "Dsuz_fpkm", "Dsuz_A", "Dmel_A", "Dbia_A", "Dtak_A", "Dsuz_X", "Dmel_X", "Dbia_X", "Dtak_X", "Dros_X", "Dros_A", "xdn", "xds", "adn", "ads", "fpkm", "Li", "Li_list", "Dsuz_A_fpkm", "Dsuz_X_fpkm", "x1m", "x2m", "x3m", "x1f", "a1m", "a1f"))
load("Dros_dnds.R")

boxplot(xdnds[,2,1], adnds[,2,1], xdnds[,4,1], adnds[,4,1], col="grey", pch=19, cex=0.8, names=c("X Dsuz-Dbia", "A Dsuz-Dbia", "X Dsuz-Dmel", "A Dsuz-Dmel"), ylim=c(0,2.5), las=1)



#Mann-Whitney U test these

#X - dN
#Dsuz-Dbia, Dsuz-Dtak
wilcox.exact(xdn[, 2, 1], xdn[, 3, 1], conf.level=0.95, conf.int=TRUE, na.rm=T)#*
#Dsuz-Dtak, Dbia-Dtak
wilcox.exact(xdn[, 3, 1], xdn[, 3, 2], conf.int=TRUE, conf.level=0.95, na.rm=TRUE)#n.s.




#X dN/dS
#Dsuz-Dbia, Dsuz-Dtak
wilcox.exact(xdnds[,2, 1], xdnds[,3,1], conf.level=0.95, conf.int=TRUE, na.rm=TRUE)
perm.test(xdnds[,2, 1], xdnds[,3,1], conf.level=0.95, conf.int=TRUE, na.rm=TRUE, exact=FALSE)
#Dsuz-Dtak, Dbia-Dtak
wilcox.exact(xdnds[,3,1], xdnds[,3,2], conf.level=0.95, conf.int=TRUE, na.rm=TRUE)#n.s

#Dsuz-Dmel, Dbia-Dmel
wilcox.exact(xdnds[,4,1], xdnds[,3,2], conf.level=0.95, conf.int=TRUE, na.rm=TRUE)#n.s
#Dsuz-Dbia, Dbia-Dmel
wilcox.exact(xdnds[,2,1], xdnds[,4,2], conf.level=0.95, conf.int=TRUE, na.rm=TRUE)#*
#Dsuz-Dbia, Dsuz-Dmel
wilcox.exact(xdnds[,2,1], xdnds[,4,1], conf.level=0.95, conf.int=TRUE, na.rm=TRUE)#*


#A -dN
#Dsuz-Dbia, Dsuz-Dtak
wilcox.exact(adn[,2, 1], adn[,3,1], conf.level=0.95, conf.int=TRUE)#*
#Dsuz-Dtak, Dbia-Dtak
wilcox.exact(adn[,3, 1], adn[,3,2], conf.level=0.95, conf.int=TRUE)#*
#Dsuz-Dmel, Dbia-Dmel
wilcox.exact(adn[,4, 1], adn[,4,2], conf.level=0.95, conf.int=TRUE)#n.s

#A dN/dS
#Dsuz-Dbia, Dsuz-Dtak
wilcox.exact(adnds[,2,1], adnds[,3,1], conf.level=0.95, conf.int=TRUE, na.rm=TRUE)#n.s
#Dsuz-Dtak, Dsuz-Dmel
wilcox.exact(adnds[,3,1], adnds[,4,1], conf.level=0.95, conf.int=TRUE)#n.s
#Dsuz-Dbia, Dsuz-Dmel
wilcox.exact(adnds[,2,1], adnds[,4,1], conf.level=0.95, conf.int=TRUE)#n.s
#Dbia-Dtak, Dbia-Dmel
wilcox.exact(adnds[,3,2], adnds[,4,2], conf.level=0.95, conf.int=TRUE)#n.s
#Dbia-Dtak, Dbia-Dmel
wilcox.exact(adnds[,3,3], adnds[,4,3], conf.level=0.95, conf.int=TRUE)#*


#where is the slow down? Grouped data show:
#Not in the autosomes, in the X.

#All data combined (MAS1, m_a, m_s)
#Dsuz-Dbia, Dsuz-Dtak
wilcox.exact(mas1[,2,1], mas1[,3,1], conf.level=0.95, conf.int=TRUE, na.rm=TRUE)#***** This is important. Just comparing overall dN/dS show a statistically significant difference between Dsuz and Dbia, but if you look at autosomal dN/dS you don't see it. The signal is in the X!!!! X slow down for Dsuz?!?!?!?!
wilcox.test(mas1[,2,1], mas1[,3,1], conf.int=TRUE, cont.level=0.95)
perm.test(mas1[,2,1], mas1[,3,1], conf.int=TRUE, conf.level=0.95,  na.rm=TRUE, exact=FALSE)

wilcox.exact(m_a[,2,1], m_a[,3,1], conf.level=0.95, conf.int=TRUE, na.rm=TRUE)


#permutation test, more conservative (I think)
#Dsuz-Dbia, Dsuz-Dtak
perm.test(adnds[,2, 1], adnds[,3,1], conf.int=TRUE, conf.level=0.95, na.rm=TRUE, exact=FALSE)#n.s
#Dsuz-Dbia, Dsuz-Dmel
perm.test(adnds[,2, 1], adnds[,4,1], conf.int=TRUE, conf.level=0.95, na.rm=TRUE, exact=FALSE)#n.s


#dN and dS for X Dsuz - Dbia, is this valid? Unequal sample sizes.




#Li's list: she requested ka/ks ratios for these genes

Li <- read.csv("~/Desktop/Molecular/SWD/SWD_genome/Li_list.csv", header=TRUE)
str(Li)
head(Li)

sum(Ortho$Gene %in% Li$Gene)


Li_list <- merge(Ortho, Li, by="Gene")
Li_list <- droplevels(Li_list)
str(Li_list)

cat(as.character(Li_list$Gene), file="Li_out.txt", sep="\r")

setwd("~/Desktop/Molecular/SWD/SWD_genome/Li_genes")


Li_files <- list.files(path="~/Desktop/Molecular/SWD/SWD_genome/Li_genes", pattern="*aligned_codon", full.names=T, recursive=FALSE)
str(Li_files)

kaks(Li_files[1])

Li_kaks <- lapply(Li_files, function(x){
	d <- read.alignment(x, format="fasta")
	out <- kaks(d)
	ka_ks <- out$ka/out$ks
	ka_ks <- as.matrix(ka_ks)
	ka_ks <- ka_ks[,c("Dsuz", "Dbia", "Dtak", "Dmel")]
	ka_ks <- ka_ks[c("Dsuz", "Dbia", "Dtak", "Dmel"),]
	return(ka_ks)
})
Li_kaks



setwd("~/Desktop/Molecular/SWD/SWD_genome/Xout")
g6851 <- read.alignment("FAMILY.6851.aligned_codon", format="fasta")
kaks(g6851)



#Breaking down Dsuz fpkm into Autosomal and X-linked genes

Dsuz_A_fpkm <- droplevels(merge(Dsuz_fpkm, Dsuz_A, by="Gene"))
str(Dsuz_A_fpkm)

Dsuz_X_fpkm <- droplevels(merge(Dsuz_fpkm, Dsuz_X, by="Gene"))
str(Dsuz_X_fpkm)

Ds_1 <- Dsuz_X_fpkm$GC3s.x[which(Dsuz_X_fpkm$GC3s.x>0.5)]

Dsuz_X_fpkm$GC3s.x[which(Dsuz_X_fpkm$GC3s.x>0.5)]


summary(Dsuz_X_fpkm$Male)#median = 9.424
male_fpkm <- Dsuz_X_fpkm$Male[which(Dsuz_X_fpkm$Male > median(Dsuz_X_fpkm$Male))]

x1m <- subset(Dsuz_X_fpkm, Male > median(Dsuz_X_fpkm$Male))
str(x1m)
boxplot(x1m$GC3s.x, Dsuz_X$GC3s, ylim=c(0,1), col="grey", names=c("High fpkm", "Normal fpkm"), las=1, pch=19, main="Male X-linked genes")

x2m <- subset(Dsuz_X_fpkm, log(Male) > 2.9)
str(x2m)
boxplot(x2m$GC3s.x, Dsuz_X$GC3s, ylim=c(0,1), col="grey", names=c("High fpkm", "Normal fpkm"), las=1, pch=19, main="Male X-linked genes")

x3m <- subset(Dsuz_X_fpkm, log(Male) > 3.5)
str(x3m)
boxplot(x3m$GC3s.x, Dsuz_X$GC3s, ylim=c(0,1), col="grey", names=c("High fpkm", "Normal fpkm"), las=1, pch=19, main="Male X-linked genes")
t.test(x3m$GC3s.x, Dsuz_X$GC3s)

x1f <- subset(Dsuz_X_fpkm, log(Female) > 3.5)
boxplot(x1f$GC3s.x, Dsuz_X$GC3s, ylim=c(0,1), col="grey", names=c("High fpkm", "Normal fpkm"), las=1, pch=19, main="Female X-linked genes")
t.test(x1f$GC3s.x, Dsuz_X$GC3s)



a1m <- subset(Dsuz_A_fpkm, log(Male) > 3.9)
str(a1m)
boxplot(a1m$GC3s.x, Dsuz_A$GC3s, ylim=c(0,1), col="grey", names=c("High fpkm", "Normal fpkm"), las=1, pch=19, main="Male autosomal genes")
t.test(a1m$GC3s.x, Dsuz_A$GC3s)
wilcox.exact(a1m$GC3s.x, Dsuz_A$GC3s, conf.level=0.95, conf.int=TRUE, na.rm=TRUE)


a1f <- subset(Dsuz_A_fpkm, log(Female) > 3.5)
str(a1f)
boxplot(a1f$GC3s.x, Dsuz_A$GC3s, ylim=c(0,1), col="grey", names=c("High fpkm", "Normal fpkm"), las=1, pch=19, main="Female autosomal genes")
t.test(a1f$GC3s.x, Dsuz_A$GC3s)



#look at male biased genes and their codon bias 
#check on X vs autosome proportion, look for underepresented X linked genes

bias <- read.csv("SWD_sex_biased.csv", header=TRUE)
str(bias)

SWD_bias <- droplevels(merge(Dsuz, bias, by="Gene"))
str(SWD_bias)

SWD_bias_m <- subset(SWD_bias, Sex_Bias = Male)
pdf(file="SWD_male_bias.pdf", bg="white")
boxplot(SWD_bias_m$GC3s, Dsuz$GC3s, col="grey", ylim=c(0,1), las=1, names=c(expression(paste(italic("D. suzukii")," male biased GC"[3])), expression(paste("Overall ", italic("D. suzukii"), " GC"[3]))), pch=19, cex=0.8)
dev.off()
wilcox.exact(SWD_bias_m$GC3s, Dsuz$GC3s, conf.level=0.95, conf.int=TRUE, na.rm=TRUE)

boxplot(SWD_bias_m$Nc, Dsuz$Nc, col="grey", las=1, names=c("Male biased GC3", "SWD X-linked GC3"), pch=19, cex=0.8)
wilcox.exact(SWD_bias_m$Nc, Dsuz$Nc, conf.level=0.95, conf.int=TRUE, na.rm=TRUE)


plot(SWD_bias_m$Fop, SWD_bias_m$GC3s)
plot(SWD_bias_m$Nc, SWD_bias_m$GC3s, pch=19, xlab="ENC", ylab=expression(paste("GC"[3])), las=1, ylim=c(0,1), xlim=c(0, 50))



#look at proportion X-linked
SWD_x_bias <- droplevels(merge(Dsuz_X, bias, by="Gene"))
str(SWD_x_bias)#65 of the biased expression genes are X-linked
SWD_a_bias <- droplevels(merge(Dsuz_A, bias, by="Gene"))
str(SWD_a_bias)#446 are autosomal

#test proportions between genome X-A and male biased X-A
x_a <- matrix(c(802, 4117, 65, 446), nrow=2)
chisq.test(x_a, simulate.p.value=TRUE, B = 10000000)
g.test(x_a)#uses Sokal and Rolhf
g.test(x_a, correct="yates")
fisher.test(x_a)

