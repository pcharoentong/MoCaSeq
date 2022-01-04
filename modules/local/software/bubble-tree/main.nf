

process bubble_tree_matched {
	tag "${meta.sampleName}"

	publishDir "${params.output_base}/${genome_build}/${meta.sampleName}/results/BubbleTree", mode: "copy"

	input:
		val (genome_build)
		val (intervals)
		tuple val (meta), path (loh_tsv), val (cn_source), path (segments_tsv)

	output:
		tuple val (meta), path ("${meta.sampleName}.Bubbletree.txt"), emit: result
		path ("*.pdf")

	when:
		!params.tiny

	script:
	"""#!/usr/bin/env Rscript

library(ggplot2)
library(BubbleTree)
library (dplyr)

# function to plot all chromosomes individually
PlotChromBT <- function(data=pred, chrom, mode="adjusted"){

  # define the data, chromosome and output name for both modes
  if(mode=="adjusted"){
    dat <- pred@rbd.adj
    plotdata <- dat[dat\$seqnames == as.character(chrom),]
    pdfname <- paste0(name,"/results/BubbleTree/",name,"_Chromosomes/adjusted/",name,".Bubbletree.adjusted.chr",chrom,".pdf")
  } else if(mode == "unadjusted"){
    dat <- pred@rbd
    plotdata <- dat[dat\$seqnames == as.character(chrom),]
    pdfname <- paste0(name,"/results/BubbleTree/",name,"_Chromosomes/unadjusted/",name,".Bubbletree.unadjusted.chr",chrom,".pdf")
  } else {
    stop("ERROR: only adjusted and unadjusted valid for mode")
  }

  # for some stupid reason, the first entry is set to color white, this way we catch it
  plotdata <- rbind(emptyElem, plotdata)

  # draw the plot
  btree <- drawBTree(btreeplotter, plotdata)

  # save the plot
  ggsave(filename = pdfname, plot = btree, width = 16, height = 9, device = "pdf")
}

add_num_probes <- function (data)
{
	if ( "Num_Probes" %in% names (data))
	{
		return (data)
	}
	else
	{
		return (data %>% dplyr::mutate (Num_Probes=(End-Start)/1000) )
	}
}

cat ("Running BubbleTree using copy number data from ${cn_source}\\n")
intervals <- strsplit ("${intervals}", ",", fixed=T)[[1]]

min_num_markers <- 10

data_loh <- read.table (file="${loh_tsv}",sep="\\t",header=T,stringsAsFactors=F)
data_cnv <- read.table (file="${segments_tsv}",sep="\\t",header=T,stringsAsFactors=F)

gr_loh <- data_loh %>%
	dplyr::mutate (start=Pos-1) %>%
	dplyr::rename (seqnames=Chrom,end=Pos,freq=Tumor_Freq) %>%
	GenomicRanges::makeGRangesFromDataFrame (keep.extra.columns=T) %>%
	GenomeInfoDb::keepSeqlevels (intervals,pruning.mode="tidy")

gr_cnv <- switch ("${cn_source}",
		"hmm-copy"={
			data_cnv %>%
			dplyr::rename (seqnames=Chrom,start=Start,end=End,seg.mean=Mean)
		},
		"dryclean-cnv-kit"={
			data_cnv %>%
			dplyr::rename (seqnames=chromosome,num.mark=probes,seg.mean=log2)
		},
		"cnv-kit-pon"={
			data_cnv
		},
		"cnv-kit-pon-cnv-kit"={
			data_cnv %>%
			dplyr::rename (seqnames=chromosome,num.mark=probes,seg.mean=log2)
		},
		{
			data_cnv
		}
	) %>%
	dplyr::select (seqnames,start,end,num.mark,seg.mean) %>%
	dplyr::filter (num.mark>=!!min_num_markers) %>%
	dplyr::mutate (seg.sig=abs(seg.mean)>0.68) %>%
	GenomicRanges::makeGRangesFromDataFrame (keep.extra.columns=T) %>%
	GenomeInfoDb::keepSeqlevels (intervals,pruning.mode="tidy")

print (gr_loh)
print (gr_cnv)

# Investigate the distribution of CNV events
cnv_plot <- ggplot (data.frame (width=width (gr_cnv),seg.mean=mcols (gr_cnv)[,"seg.mean"],seg.sig=mcols (gr_cnv)[,"seg.sig"]),aes (seg.mean, width)) +
	scale_y_continuous (trans="log10") +
	scale_x_continuous (limits=c(-5,5),oob=scales::squish) +
	geom_point (aes(colour=seg.sig)) +
	theme_bw ()

ggsave(filename="${meta.sampleName}.Bubbletree.cnv.pdf", plot = cnv_plot, width = 16, height = 9, device = "pdf")

# Get rid of the unreliable segments
mcols(gr_cnv)[!mcols(gr_cnv)[,"seg.sig"],"seg.mean"] <- 0

# for some stupid reason, the first entry is set to color white, this way we catch it
emptyElem <- data.frame(seqnames=0,start=0, end=0,width=0, strand="*", seg.id=1000, num.mark=1000, lrr=0, kurtosis=0, hds=0, hds.sd=0,het.cnt=0,seg.size=0)

r <- new("RBD", unimodal.kurtosis=-0.1)

rbd <- makeRBD(r, gr_loh,gr_cnv)
btreepredictor <- new("BTreePredictor", rbd=rbd, max.ploidy=6, prev.grid=seq(0.2,3, by=0.01))
btreepredictor@config\$min.segSize <- ifelse(max(btreepredictor@rbd\$seg.size,na.rm=TRUE) < 0.4, 0.1, 0.4)
# ^was copied from the manual, however the vignette filters on num.markers instead
# rbd\$seg.size <- rbd\$num.mark / total.mark * 100
# The default: min.segSize = 0.5
# for the PON CNVKit calls the top call has seg.size ~ 1.5, therefore the top 65 calls would be used
# There is a problem if there are only a handful of calls above the 0.4 threshold
btreepredictor@config\$min.segSize <- ifelse(max(btreepredictor@rbd\$seg.size,na.rm=TRUE) < 0.4 || nrow (btreepredictor@rbd[btreepredictor@rbd[,"seg.size"] > 0.4,]) < 50, 0.1, 0.4)
cat ("\\nUsing min.segSize: ",btreepredictor@config\$min.segSize," retained ",nrow (btreepredictor@rbd[btreepredictor@rbd[,"seg.size"] > btreepredictor@config\$min.segSize,])," segments\\n")
pred <- btpredict(btreepredictor)

info <- info(pred)
write.table(info,"${meta.sampleName}.Bubbletree.txt", col.names=F, row.names=F, quote=F, sep="\\t")

cat("\\nPurity/Ploidy: ", info, "\\n")

# plot full unadjusted
btreeplotter <- new("BTreePlotter", branch.col="gray50")
pred@rbd <- rbind(emptyElem, pred@rbd) # for some stupid reason, the first entry is set to color white, this way we catch it
btree <- drawBTree(btreeplotter, pred@rbd)
ggsave(filename="${meta.sampleName}.Bubbletree.unadjusted.pdf", plot = btree, width = 16, height = 9, device = "pdf")

# plot full adjusted
btreeplotter <- new("BTreePlotter", branch.col="gray50")
pred@rbd.adj <- rbind(emptyElem, pred@rbd.adj) # for some stupid reason, the first entry is set to color white, this way we catch it
btree <- drawBTree(btreeplotter, pred@rbd.adj)
ggsave(filename="${meta.sampleName}.Bubbletree.adjusted.pdf", plot = btree, width = 16, height = 9, device = "pdf")

# plot shift from unadjusted to adjusted
btreeplotter <- new("BTreePlotter", branch.col="gray50")
rbd1 <- pred@rbd
rbd2 <- pred@rbd.adj
arrows <- trackBTree(btreeplotter, rbd1, rbd2, min.srcSize=0.01, min.trtSize=0.01)
btree <- drawBTree(btreeplotter, rbd2) + arrows
ggsave(filename="${meta.sampleName}.Bubbletree.adjustmentshift.pdf", plot = btree, width = 16, height = 9, device = "pdf")


	"""

	stub:
	"""#!/usr/bin/env bash

if [[ "${params.stub_json_map?.bubble_tree_matched}" == "null" ]]; then
	cp ${params.stub_dir}/${genome_build}/${meta.sampleName}/results/BubbleTree/${meta.sampleName}.Bubbletree.txt .
fi

touch ${meta.sampleName}.Bubbletree.txt
touch ${meta.sampleName}.Bubbletree.adjusted.pdf
touch ${meta.sampleName}.Bubbletree.adjustmentshift.pdf
touch ${meta.sampleName}.Bubbletree.cnv.pdf
touch ${meta.sampleName}.Bubbletree.unadjusted.pdf
	"""

}



