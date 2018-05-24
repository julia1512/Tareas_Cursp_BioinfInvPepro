rm(list=ls())
#Definir region
chrom = "chr6"
chromstart = 31541312
chromend = 31541551

#Plot de Genes anotados
library(Sushi)
Genes.bed<-read.delim("fig/gene_lta_gene_2.txt",header = F)
head(Genes.bed)
par(mar=c(4,3.5,2,3))
plotGenes(Genes.bed,chrom,chromstart,chromend ,
          types=Genes.bed$types,
          maxrows=1,bheight=0.2,
          plotgenetype="box",bentline=FALSE, #plotgenetype="arrow",bentline=FALSE,
          labeloffset=.4,fontsize=1.2,arrowlength = 0.025,
          labeltext=F, rowlabels = "B.", rowlabelcol = "black",col = "DarkBlue",packrow = T)
labelgenome(chrom,chromstart,chromend,n=7,scale="Mb")

#Plot de metilación en CpGs
meth_huvec.bedgraph<-read.delim("fig/methylation_huvec.txt",he=F)
meth_GM.bedgraph<-read.delim("fig/methylation_GM1287.txt",he=F)
meth_K562c.bedgraph<-read.delim("fig/methylation_K562.txt",he=F)
meth.bedgraph<-rbind(meth_GM.bedgraph,meth_K562c.bedgraph)
meth.bedgraph<-rbind(meth.bedgraph,meth_huvec.bedgraph)
row_cell<-c(1,1,2,2,3,3)
meth.bedgraph<-cbind(meth.bedgraph,row_cell)
par(mar=c(4,6,2,3))
plotBed(beddata = meth.bedgraph,chrom = chrom,
        chromstart = chromstart,chromend =chromend,
        rownumber = meth.bedgraph$row_cell, type = "region", color=c("#7F0F7E","#FE7F23","#7F0F7E","#FE7F23","#FE7F23","#FE7F23"),
        row="given", plotbg="white",
        rowlabels=unique(c("GM1287","K562","HUVEC")), 
        rowlabelcol="black",
        rowlabelcex=1)
labelgenome(chrom,chromstart,chromend,n=7,scale="Mb")

