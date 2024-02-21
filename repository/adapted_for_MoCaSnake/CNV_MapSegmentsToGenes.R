#!/usr/bin/Rscript

##########################################################################################
##
## CNV_MapSegmentsToGenes.R
##
## Takes the segment file from either HMMCopy or Copywriter and maps them to genes.
## Customized for MoCaSnake such that input and output directory can be provided by Snakemake
##
##########################################################################################

args <- commandArgs(TRUE)

segment_file = args[1]
name = args[2]
genecode_file = args[3]
resolution = args[4]
CGC = args[5]
TruSight = args[6]
output_folder = args[7]

suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(data.table))

genesDT <- readRDS(genecode_file)
genesDT <- data.table(genesDT)
genesGR <- makeGRangesFromDataFrame(genesDT, keep.extra.columns = T)

AnnotateSegment <- function(segDF){
  segDF[,"chr"]=paste0("chr",segDF[,"chr"])
  segGR <- makeGRangesFromDataFrame(segDF, keep.extra.columns = T)
  hits <- findOverlaps(segGR, genesGR)
  returnDat <- genesDT[subjectHits(hits), .(chr, start, end, geneName, geneID)] # more stuff
  return(data.frame(returnDat))
}

cnv <- data.frame(Name=NULL,Chrom=NULL, Start=NULL, End=NULL, Mean=NULL,Gene=NULL)
segment <- fread(segment_file)

segment[Chrom == 23, Chrom := "X"]
segment[Chrom == 24, Chrom := "Y"]
genesDT[chr == 23, Chrom := "X"]
genesDT[chr == 24, Chrom := "Y"]
chromorder <- c(1:22, "X", "Y")
segment <- data.frame(segment)

for (i in 1:nrow(segment)) {
temp=NULL
temp=as.data.frame(segment[i,c("Chrom","Start","End")])
colnames(temp)=c("chr","start","end")
results=AnnotateSegment(temp)
    if (nrow(results) > 0) 
    {
        colnames(results)=c("Chrom", "Start", "End", "Gene", "GeneID")
        results[,"Chrom"]=gsub("chr","",results[,"Chrom"])
        results$Mean=segment[i,"Mean"]
        results$Name=name
        results=results[,c("Name", "Chrom", "Start", "End", "Mean","Gene", "GeneID")]
        cnv=rbind(cnv,results)
    }
}

cnv=tbl_df(cnv) %>% mutate_each(as.character)

cnv = cnv %>% 
  group_by(Gene) %>% 
  arrange(desc(abs(as.numeric(Mean))),.by_group=T) %>% 
  filter(row_number()==1) %>% 
  arrange(as.numeric(Start),as.numeric(End))

# sort by chromosome
cnv$Chrom <- factor(cnv$Chrom, levels = chromorder)
cnv$Start <- as.numeric(cnv$Start)
setorder(cnv, Chrom, Start)

write.table(cnv,paste(output_folder,"/",name,".",resolution,".genes.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

CGC=read.delim(CGC,header=T,sep="\t")

cnv_cgc = cnv %>%
  filter(Gene %in% CGC[,1])

write.table(cnv_cgc,paste(output_folder,"/",name,".",resolution,".genes.CGC.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

cnv_cgc = cnv_cgc %>%
  filter(abs(as.numeric(Mean)) > 0.75)

write.table(cnv_cgc,paste(output_folder,"/",name,".",resolution,".genes.OnlyImpact.CGC.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

TruSight=read.delim(TruSight,header=T,sep="\t")

cnv_ts = cnv %>%
filter(Gene %in% TruSight[,1])

write.table(cnv_ts,paste(output_folder,"/",name,".",resolution,".genes.TruSight.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

cnv_ts = cnv_ts %>%
filter(abs(as.numeric(Mean)) > 0.75)

write.table(cnv_ts,paste(output_folder,"/",name,".",resolution,".genes.OnlyImpact.TruSight.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

