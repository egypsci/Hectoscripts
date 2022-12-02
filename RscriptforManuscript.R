## This script uses the gene expression matrices matrix_uA and matrix_mA.
## Ahmed Elewa elewa@ub.edu
## 11/13/2022
##  
## 
### Script for Subramanian et al. Dev Cell ###

source('WABI/Rscripts/functions.R')
source('WABI/Rscripts/pca_plotting_v2.R')

library(plotrix)

# Unzip the following files first
df_uA <- read.table('matrix_uA',header=T) # uniquely mapped reads
df_mA <- read.table('matrix_mA',header=T) # multimapped reads
df_aA <- df_uA + df_mA

facs <- read.csv('facsdata.csv',row.names = 1,header=TRUE)

# Identify data for ERCC spike-ins and non-coding RNAs
rna_aA  <-grep("Infernal",rownames(df_aA))
ercc_aA <-grep("^ERCC-[[:digit:]]",rownames(df_aA))

# Number of detected genes and ERCC spike-ins
nD.genes      <- apply(df_aA[-c(ercc_aA,rna_aA),],2,detect.genes)
nD.ercc       <- apply(df_aA[ercc_aA,],2,detect.genes)

# Filter out cells with less than 10 ERCC
passed <- (nD.ercc >=10)

# Colors
popAcolor <- rgb((252/255),(169/255),(133/255),0.5)
popBcolor <- rgb((133/255),(202/255),(93/255),0.5)
ribocolor <- rgb((1/255),(100/255),(255/255),0.5)
TEcolor   <- rgb((255/255),(1/255),(100/255),0.5)

# Keep protein-coding genes with more than one mapped read
df_aA2 <- df_aA[-c(ercc_aA,rna_aA),]
df_aA2 <- df_aA2[rownames(df_aA2)!='UNKNOWN',]
df_aA2 <- df_aA2[rownames(df_aA2)!='unknown',]
keep<-which(rowSums(df_aA2)>1)

# Normalize using kept genes and passed cells
df_aAn <- sweep(df_aA2[keep,passed], 2, colSums(df_aA2[keep,passed]), FUN="/") *1000000

# Define two populations
popA <- colnames(df_aA2[keep,nD.genes >= 2000 & nD.ercc >=10])
popB <- colnames(df_aA2[keep,nD.genes < 2000  & nD.ercc >=10])

# Genes of interest
ribogenes <- c('RL10','RL10A','RL10L','RL11','RL12','RL13','RL13A','RL14','RL15','RL17','RL18','RL18A',
               'RL19','RL21','RL22','RL22L','RL23','RL23A','RL24','RL26','RL26L','RL27','RL27A','RL28',
               'RL29','RL3','RL30','RL31','RL32','RL34','RL35','RL35A','RL36','RL36A','RL37','RL37A','RL38',
               'RL39','RL4','RL40','RL4B','RL5','RL5B','RL6','RL7','RL7A','RL7L','RL8','RL9','RLA0','RLA1',
               'RLA2','RS10','RS11','RS12','RS13','RS14','RS15','RS15A','RS16','RS17','RS18','RS19','RS2',
               'RS20','RS21','RS23','RS24','RS25','RS26','RS27','RS27A','RS27L','RS28','RS29','RS3','RS30',
               'RS31','RS3A','RS3AB','RS4','RS4X','RS4Y1','RS5','RS6','RS7','RS8','RS9')
TEs <- c('DPOL','ENR1','ENV','ERVV1','ERVV2','GAG','HARB1','LIN1','LITD1','LORF1','LORF2','POL',
         'POL1','POL4','RTL1','Transposon_cacta','Transposon_Crypton','Transposon_DDE_1','Transposon_gypsy',
         'Transposon_hAT','Transposon_helitronORF','Transposon_ISb','Transposon_ISC1316','Transposon_LINE',
         'Transposon_ltr_Roo','Transposon_mariner','Transposon_MuDR_A_B','Transposon_P_element',
         'Transposon_piggybac','Transposon_TY1_Copia')


# PCA (all cells)
PC_aAn <-run.pca(df_aAn)

# PC contributions and top contributers to each PC (all cells)
PC_aAn.vars<- PC_aAn$sdev^2
PC_aAn.vars<- PC_aAn.vars/sum(PC_aAn.vars)
PC_aAn.aload <- abs(PC_aAn$rotation[,1:10])
PC_aAn.contr<- sweep(PC_aAn.aload, 2, colSums(PC_aAn.aload), "/")

# PCA (popB)
PC_aAnB <-run.pca(df_aAn[,popB])

# PC contributions and top contributers to each PC (popB)
PC_aAnB.vars<- PC_aAnB$sdev^2
PC_aAnB.vars<- PC_aAnB.vars/sum(PC_aAnB.vars)
PC_aAnB.aload <- abs(PC_aAnB$rotation[,1:10])
PC_aAnB.contr<- sweep(PC_aAnB.aload, 2, colSums(PC_aAnB.aload), "/")

# PCA (popA)
PC_aAnA <-run.pca(df_aAn[,popA])

# PC contributions and top contributers to each PC (popB)
PC_aAnA.vars<- PC_aAnA$sdev^2
PC_aAnA.vars<- PC_aAnA.vars/sum(PC_aAnA.vars)
PC_aAnA.aload <- abs(PC_aAnA$rotation[,1:10])
PC_aAnA.contr<- sweep(PC_aAnA.aload, 2, colSums(PC_aAnA.aload), "/")

# Negative correlation between ribogenes and TEs
cor.test(colMeans(df_aAn[TEs,]),colMeans(df_aAn[ribogenes,]))
cor.test(colMeans(df_aAn[TEs,popA]),colMeans(df_aAn[ribogenes,popA]))
cor.test(colMeans(df_aAn[TEs,popB]),colMeans(df_aAn[ribogenes,popB]))

cor.test(colMeans(df_aAn[TEs,popB][substr(popB,1,2)=='P1']),colMeans(df_aAn[ribogenes,popB][substr(popB,1,2)=='P1']))
substr(popB,1,2)

############################
##        Figures         ##
############################


# Fig 1b.
plot(nD.ercc[passed],nD.genes[passed],main='',
     pch=20,
     col=ifelse(names(nD.genes[passed])%in%popA,rgb((252/255),(169/255),(133/255),0.5),rgb((133/255),(202/255),(93/255),0.5)),
     ylim=c(0,8500),ylab = 'detected genes',xlab='detected ERCC', las=1,cex.axis=1)
legend('topleft',legend=c('kiloscript','hectoscript'),
       col=(c(rgb((252/255),(169/255),(133/255),0.5),rgb((133/255),(202/255),(93/255),0.5))),pch=20,bty='n')


# Fig 1c. 
plot(PC_aAn$x[,'PC1'],nD.genes[passed],pch=20,
     col=ifelse(names(nD.genes[passed])%in%popA,popAcolor,popBcolor),
     main='',las=2,cex.axis=1,ylab='detected genes',xlab='PC1')
abline(v=0,lty=3)
legend('topleft',legend='kiloscript',bty='n',cex=1.4)
legend('topright',legend='hectoscript',bty='n',cex=1.4)

# Fig 1f. MKNK2 is a marker for popB
plot(PC_aAn$x[,'PC1'],as.numeric(df_aAn['MKNK2',]),
     col=ifelse(names(nD.genes[passed])%in%popA,popAcolor,popBcolor),yaxt='n',
     pch=20,main='',cex.axis=0.7,ylab='MKNK2 expression CPM/1000',xlab='PC1')
axis(2,at=c(50000,100000,150000,200000),labels=c('50','100','150','200'),las=2,cex.axis=0.7)
abline(v=0,lty=3)
text(-200,200000,'Kiloscript')
text(105,200000,'Hectoscript')

# Fig 1d.
nPlot <- 20 
for (i in 1:1) {
  top<-order(PC_aAn$rotation[,i],decreasing=T)[1:nPlot]
  bottom<-order(PC_aAn$rotation[,i],decreasing=F)[1:nPlot]
  barplot(PC_aAn.contr[top,i],main='',ylab='',
          las=2,horiz=T,cex.axis=0.5,cex.names=0.5)
}

# Fig 1e.
boxplot(as.numeric(df_aAn['SMUF2',]),as.numeric(df_aAn['MKNK2',]),as.numeric(df_aAn['TPD52',]),
        as.numeric(df_aAn['WDR26',]),as.numeric(df_aAn['QRIC1',]),as.numeric(df_aAn['CHSS1',]),
        as.numeric(df_aAn['MARCS',]),as.numeric(df_aAn['KS6B1',]),outline = F,las=2,
        ylab='Hectoscript expression CPM',xaxt='n',ylim=c(0,200000),frame.plot=T,cex.axis=0.7)
axis(1,at=1:8,labels = c('SMURF2','MKNK2','TPD52','WDR26','QRIC1','CHSS1','MARCS','KS6B1'),las=2,cex.axis=0.7)



# Fig 1g.
pca.plot(PC_aAnB,pch=16,main='',
         col=color.scale(log2(colMeans(df_aAn[TEs,popB])),c(0,1,1),c(1,1,0),0),xrange(0,14),las=2,cex.axis=1)

pca.plot(PC_aAnB,pch=16,main='',
         col=color.scale(log(colMeans(df_aAn[ribogenes,popB])),c(0,1,1),c(1,1,0),0),xrange(0,14),las=2,cex.axis=1)

# legend
par(mar = c(10,10,14,10))
plot(0:14,rep(1,15),pch=15,ylim=c(1,100),cex=10,col=color.scale(0:14,c(0,1,1),c(1,1,0),0),yaxt='n',bty='n',ylab='',xlab='')


# Fig 1h.
plot(PC_aAn$x[,'PC1'],log(colMeans(df_aAn[ribogenes,])),
     ylim=c(4.5,10.5),
     col=ribocolor,
     pch=20,main='',las=1,cex.axis=0.7,ylab='ln(meanCPM)',xlab='PC1')

points(PC_aAn$x[,'PC1'],log(colMeans(df_aAn[TEs,])),
       col=TEcolor,pch=16)
splineTE <- smooth.spline(PC_aAn$x[,'PC1'],log(colMeans(df_aAn[TEs,])), df=4)
predictTE <- predict(splineTE, PC_aAn$x[,'PC1'])
splineRibogene <- smooth.spline(PC_aAn$x[,'PC1'],log(colMeans(df_aAn[ribogenes,])), df=4)
predictRibogene <- predict(splineRibogene, PC_aAn$x[,'PC1'])
lines(splineTE,col=TEcolor,lwd=6)
lines(splineRibogene,col=ribocolor,lwd=6)
legend('topleft',pch=c(20,16),col=c(ribocolor,TEcolor),c('Ribogene','TE'),bty='n')
abline(v=0,lty=3)



## Supplemetnary Figure
facs2 <- facs[rownames(facs)%in%popA|rownames(facs)%in%popB, ]

colbatch = c('chartreuse3','darkgoldenrod')
pca.plot(PC_aAn,pch=16,main='PCA and Batch Effect',col=colbatch[as.numeric(substr((colnames(df_aAn)),2,2))],
         ylim=c(-220,220),xlim=c(-420,300),las=1)


nPlot <- 20 #number of genes top plot per pc.
for (i in 1:1) {
  top<-order(PC_aAn$rotation[,i],decreasing=T)[1:nPlot]
  bottom<-order(PC_aAn$rotation[,i],decreasing=F)[1:nPlot]
  #barplot(PC_aAn.contr[top,i],main='',ylab='',
  #        las=2,horiz=T,cex.axis=0.5,cex.names=0.5)
  barplot(PC_aAn.contr[bottom,i],main='Top Contributors PC1-ve',ylab='',
          las=2,horiz=T,cex.axis=1,cex.names=1)
}


par(mfrow=c(2,2))
plot(nD.ercc[passed],nD.genes[passed],main='',
     pch=20,
     col=ifelse(names(nD.genes[passed])%in%popA,rgb((252/255),(169/255),(133/255),0.5),rgb((133/255),(202/255),(93/255),0.5)),
     ylim=c(0,8500),ylab = 'detected genes',xlab='detected ERCC', las=1,cex.axis=1)
legend('topleft',legend=c('kiloscript','hectoscript'),
       col=(c(rgb((252/255),(169/255),(133/255),0.5),rgb((133/255),(202/255),(93/255),0.5))),pch=20,bty='n')

plot(facs2$FSC,facs2$SSC,pch=20,
     col=ifelse(rownames(facs2)%in%popA,rgb((252/255),(169/255),(133/255),0.5),rgb((133/255),(202/255),(93/255),0.5)),
     xlab='FSC (size)',ylab='SSC (internal complexity)',las=1)
boxplot(facs2[rownames(facs2)%in%popA,'FSC'],facs2[rownames(facs2)%in%popB,'FSC'],
        names=c('kiloscripts','hectoscripts'),main='FSC (size)')
boxplot(facs2[rownames(facs2)%in%popA,'SSC'],facs2[rownames(facs2)%in%popB,'SSC'],
        names=c('kiloscripts','hectoscripts'),main='SSC (internal complexity)')

cor.test(facs2[rownames(facs2)%in%popA,'FSC'],head(facs2[rownames(facs2)%in%popB,'FSC'][order(facs2[rownames(facs2)%in%popB,'FSC'],decreasing=FALSE)],350))
cor.test(facs2[rownames(facs2)%in%popA,'SSC'],head(facs2[rownames(facs2)%in%popB,'SSC'][order(facs2[rownames(facs2)%in%popB,'SSC'],decreasing=FALSE)],350))

t.test(facs2[rownames(facs2)%in%popA,'FSC'],facs2[rownames(facs2)%in%popB,'FSC'])
t.test(facs2[rownames(facs2)%in%popA,'SSC'],facs2[rownames(facs2)%in%popB,'SSC'])

