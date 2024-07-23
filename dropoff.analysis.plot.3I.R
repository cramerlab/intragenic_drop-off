########################################################################################
# FACT maintains chromatin architecture and thereby stimulates RNA polymerase II pausing 
# during transcription in vivo

# Kristina Žumer et al., Mol Cell, 2024
# DOI:https://doi.org/10.1016/j.molcel.2024.05.003

# Script for estimation of drop-off probability from TT-seq data
# Author : Arjun Devadas, MPI-NAT
########################################################################################

# load required packages
{
  library(rtracklayer)
  library(GenomicAlignments)
  library(Rsamtools)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggpubr)
  library(magrittr)
  library(doParallel)
}

# paths and input variables
{
  setwd("") # change accordingly
  
  # anno
  anno.path = "data/FACT_expressed_non_ovelapping_gene_annotation.RData"
  anno = get(load(anno.path)) # granges object
  
  # cov_i and cov_f annotation
  start.width = 3e3
  anno = anno[which(width(anno) > (start.width + 5e2))]
  coverage.initial.anno = promoters(anno, upstream = 0, downstream = start.width) # beginning of genes annotation
  coverage.final.anno = resize(anno, fix = "end", width = (width(anno) - start.width)) # gene body annotation
  
  # bam files
  bam.files.path = ""
  bam.files = list.files(bam.files.path, pattern = ".bam$", full.names = F) # 8 samples, 4 conditions 2 replicates each
  
 # chromosome names and lengths
  chrs.lengths.path = "data/refseq_chrs_lengths.RData"
  chrs.lengths = get(load(chrs.lengths.path))
  
  # antisense bias ratio
  antisense.bias.ratio = "data/FACT_antisense_bias_ratio.RData"
  
  # outputs
  out.file.path = ""
}

# functions
{
  # create.transcribed.bases.rle.tracks.anno
  create.transcribed.bases.rle.tracks = function(bam.files,
                                                 bam.input.folder,
                                                 rle.out.folder,
                                                 prefix,
                                                 human.chrs,
                                                 human.chrs.lengths,
                                                 strand.specific = TRUE,
                                                 remove.duplicates = FALSE,
                                                 size.selection = FALSE,
                                                 strand.mode = 1)
  {
    dir.create(file.path(rle.out.folder, "transcribed.bases.rle.tracks"))
    for (bam.file in bam.files){
      dir.create(file.path(file.path(rle.out.folder, "transcribed.bases.rle.tracks", bam.file)))
      print(bam.file)
      registerDoParallel(cores = mc.cores)
      build.transcribed.bases.coverage.track.list = function(which.chr){
        transcribed.bases.coverage.track.list.chr = list()
        param = ScanBamParam(which=GRanges(seqnames = which.chr,ranges = IRanges(0,human.chrs.lengths[which.chr])))
        bam = readGAlignmentPairs(file = file.path(paste0(bam.input.folder, bam.file)),param = param, strandMode = strand.mode)
        bam = bam[start(left(bam)) <= end(right(bam))]
        bam = bam[which(seqnames(bam) == which.chr)] # Remove non-concordant mappings
        if(length(bam) != 0){
          if(size.selection) {
            bam = bam[(end(right(bam)) - start(left(bam))) <= 500] # Size selection
          }
          if(remove.duplicates) {
            bam = bam[!duplicated(paste(start(left(bam)),end(right(bam)))),] # Remove duplicates
          }
          starts = end(left(bam[strand(bam) == "+"]))
          transcribed.bases.coverage.track.list.chr[["+"]] = 0
          if(length(starts) != 0){
            rle.vec = Rle(0,human.chrs.lengths[which.chr])
            coverage.vec = coverage(GRanges(seqnames = which.chr,ranges = IRanges(start = start(left(bam[strand(bam) == "+"])),end = end(right(bam[strand(bam) == "+"])))))[[which.chr]]
            rle.vec[1:length(coverage.vec)] = coverage.vec
            transcribed.bases.coverage.track.list.chr[["+"]] = rle.vec
          }

          starts = end(left(bam[strand(bam) == "-"]))
          transcribed.bases.coverage.track.list.chr[["-"]] = 0
          if(length(starts) != 0){
            rle.vec = Rle(0,human.chrs.lengths[which.chr])
            coverage.vec = coverage(GRanges(seqnames = which.chr,ranges = IRanges(start = start(left(bam[strand(bam) == "-"])),end = end(right(bam[strand(bam) == "-"])))))[[which.chr]]
            rle.vec[1:length(coverage.vec)] = coverage.vec
            transcribed.bases.coverage.track.list.chr[["-"]] = rle.vec
          }
          if(!strand.specific) {
            transcribed.bases.coverage.track.list.chr[["+"]] = transcribed.bases.coverage.track.list.chr[["+"]] + transcribed.bases.coverage.track.list.chr[["-"]]
            transcribed.bases.coverage.track.list.chr[["-"]] = transcribed.bases.coverage.track.list.chr[["+"]]
          }


          save(transcribed.bases.coverage.track.list.chr,file = file.path(rle.out.folder, "transcribed.bases.rle.tracks", bam.file, paste0(prefix, which.chr, ".RData")))
          return()
        }else{
          transcribed.bases.coverage.track.list.chr[["+"]] = 0
          if(!strand.specific){
            transcribed.bases.coverage.track.list.chr[["+"]] = 0
            transcribed.bases.coverage.track.list.chr[["-"]] = 0
          }
          save(transcribed.bases.coverage.track.list.chr,file = file.path(rle.out.folder, "transcribed.bases.rle.tracks", bam.file, paste0(prefix, which.chr, ".RData")))
          return()
        }

      }

      transcribed.bases.coverage.track.list = foreach(n = human.chrs,.noexport = setdiff(ls(),c("human.chrs.lengths"))) %dopar% build.transcribed.bases.coverage.track.list(n)
    }
  }

  # calculate.counts
  calculate.counts = function(anno,
                              rle.location,
                              bam.files,
                              rle.prefix,
                              out.folder,
                              file.name,
                              antisense.bias.ratio.location = "",
                              chrs,
                              chrs.lengths,
                              antisense = F,
                              find.antisense = T)
  {
    if (TRUE){
      counts = list()
      for (bam.file in bam.files){
        print(bam.file)
        index.subsets = split(1:nrow(anno),paste(as.character(anno[,"strand"]),"/",anno[,"chr"],sep = ""))
        coverage.list = list()

        build.fragment.counts.list = function(j){
          from.transcript = strand.chr.anno[j,"start"]
          to.transcript = strand.chr.anno[j,"end"]
          sum.transcript = sum(as.vector(strand.chr.fragment.counts.from.bam[from.transcript:to.transcript]))
          names(sum.transcript) = j
          return(sum.transcript)
        }

        for (index.subset in names(index.subsets)){
          print(index.subset)
          fragment.mid.track.list.chr = get(load(file.path(paste0(rle.location,bam.file,"/",rle.prefix,unlist(strsplit(index.subset,split = "/"))[2],".RData"))))

          if(unlist(strsplit(index.subset,split = '/'))[1] == "*"){
            strand.chr.fragment.counts.from.bam = (fragment.mid.track.list.chr[["+"]] + fragment.mid.track.list.chr[["-"]])/2
          }
          else{
            strand.chr.fragment.counts.from.bam = fragment.mid.track.list.chr[[unlist(strsplit(index.subset,split = '/'))[1]]]
          }

          strand.chr.anno = anno[index.subsets[[index.subset]],c("start","end","width")]

          registerDoParallel(cores = mc.cores)
          coverage.list = c(coverage.list, foreach(n = rownames(strand.chr.anno),.noexport = setdiff(ls(),c("strand.chr.fragment.counts.from.bam","strand.chr.anno","build.fragment.counts.list"))) %dopar% build.fragment.counts.list(n))
        }
        counts[[bam.file]] = unlist(coverage.list,recursive = TRUE,use.names = TRUE)
      }
      counts = sapply(counts,c)
      print("counts created")
      #rownames(counts) = as.character(anno[as.numeric(unlist(index.subsets,recursive = TRUE,use.names = TRUE)),"trid"])
      colnames(counts) = bam.files
      if(antisense) {
        # Antisense bias
        if(find.antisense){
          antisense.bias.ratio.mat = cbind()
          for (bam.file in bam.files){
            registerDoParallel(cores = mc.cores)
            build.antisense.bias.ratios = function(which.chr,distance = 0){
              chr.coverage.list = get(load(file.path(paste0(rle.location,bam.file,"/",rle.prefix,which.chr,".RData"))))
              anno.positions = Rle(0,chrs.lengths[which.chr])
              check.ids = rownames(anno[which(anno[,"strand"] == "+" & anno[,"chr"] == which.chr),])
              for (check.id in check.ids){anno.positions[(anno[check.id,"start"] - distance):(anno[check.id,"end"] + distance)] = 1}
              diff.positions = Rle(1,chrs.lengths[which.chr])
              check.ids = rownames(anno[which(anno[,"strand"] == "-" & anno[,"chr"] == which.chr),])
              for (check.id in check.ids){diff.positions[(anno[check.id,"start"] - distance):(anno[check.id,"end"] + distance)] = 0}
              positions = anno.positions*diff.positions
              runValue(positions) = as.logical(runValue(positions))
              sense = chr.coverage.list[["+"]][positions]
              antisense = chr.coverage.list[["-"]][positions]
              anno.positions = Rle(0,chrs.lengths[which.chr])
              check.ids = rownames(anno[which(anno[,"strand"] == "-" & anno[,"chr"] == which.chr),])
              for (check.id in check.ids){anno.positions[(anno[check.id,"start"] - distance):(anno[check.id,"end"] + distance)] = 1}
              diff.positions = Rle(1,chrs.lengths[which.chr])
              check.ids = rownames(anno[which(anno[,"strand"] == "+" & anno[,"chr"] == which.chr),])
              for (check.id in check.ids){diff.positions[(anno[check.id,"start"] - distance):(anno[check.id,"end"] + distance)] = 0}
              positions = anno.positions*diff.positions
              runValue(positions) = as.logical(runValue(positions))
              sense = c(sense,chr.coverage.list[["-"]][positions])
              antisense = c(antisense,chr.coverage.list[["+"]][positions])
              sense[which(sense < 100)] = NA
              antisense[which(antisense == 0)] = NA
              return(median(antisense/sense,na.rm = TRUE))
            }
            antisense.bias.ratios = as.vector(foreach(n = chrs) %dopar% build.antisense.bias.ratios(n))

            antisense.bias.ratio.mat = cbind(antisense.bias.ratio.mat,antisense.bias.ratios)
          }

          antisense.bias.ratio.mat = apply(antisense.bias.ratio.mat,c(1,2),as.numeric)
          rownames(antisense.bias.ratio.mat) = chrs
          colnames(antisense.bias.ratio.mat) = bam.files

          antisense.bias.ratio = apply(antisense.bias.ratio.mat,2,function(x){median(x,na.rm = TRUE)})
          antisense.bias.ratio[is.na(antisense.bias.ratio)] = 0
          save(antisense.bias.ratio, file = antisense.bias.ratio.location)
        }
        # Antisense counts
        antisense.counts = list()
        for (bam.file in bam.files) {
          index.subsets = split(1:nrow(anno),paste(Vectorize(strand.switch)(as.character(anno[,"strand"])),"/",anno[,"chr"],sep = ""))
          coverage.list = list()

          build.fragment.counts.list = function(j){
            from.transcript = strand.chr.anno[j,"start"]
            to.transcript = strand.chr.anno[j,"end"]
            sum.transcript = sum(as.vector(strand.chr.fragment.counts.from.bam[from.transcript:to.transcript]))
            names(sum.transcript) = j
            return(sum.transcript)
          }

          for (index.subset in names(index.subsets)){
            fragment.mid.track.list.chr = get(load(file.path(paste0(rle.location,bam.file,"/",rle.prefix,unlist(strsplit(index.subset,split = "/"))[2],".RData"))))

            strand.chr.fragment.counts.from.bam = fragment.mid.track.list.chr[[unlist(strsplit(index.subset,split = "/"))[1]]]
            strand.chr.anno = anno[index.subsets[[index.subset]],c("start","end","width")]

            registerDoParallel(cores = mc.cores)
            coverage.list = c(coverage.list, foreach(n = rownames(strand.chr.anno),.noexport = setdiff(ls(),c("strand.chr.fragment.counts.from.bam","strand.chr.anno","build.fragment.counts.list"))) %dopar% build.fragment.counts.list(n))
          }
          antisense.counts[[bam.file]] = unlist(coverage.list,recursive = TRUE,use.names = TRUE)
        }
        antisense.counts = sapply(antisense.counts,c)
        print("antisense counts created")
        #rownames(antisense.counts) = as.character(anno[as.numeric(unlist(index.subsets,recursive = TRUE,use.names = TRUE)),"trid"])
        colnames(antisense.counts) = bam.files

        antisense.bias.ratio = get(load(antisense.bias.ratio.location))
        antisense.bias.ratio = antisense.bias.ratio[bam.files]

        if(all(colnames(counts) == names(antisense.bias.ratio))) {
          print("saving counts")
          index = rownames(anno)
          counts.antisense.corrected = t(t(counts[index, ] - t(t(antisense.counts[index, ])*antisense.bias.ratio))/(1 - antisense.bias.ratio^2))
          counts.antisense.corrected[counts.antisense.corrected < 0] = 0
          counts = cbind(counts[index, ], antisense.counts[index, ])
          save(counts, file = paste0(out.folder, file.name,".RData"))
          save(counts.antisense.corrected, file = paste0(out.folder,file.name, ".antisense.corrected.RData"))
        }
        else {
          print("Incorrect antisense bias ratio file!!")
        }
      }
      else {
        print("saving counts")
        index = rownames(anno)
        counts = cbind(counts[index, ])
        save(counts, file = paste0(out.folder, file.name,".RData"))
      }
    }
  }

  # additional functions
  right = function(x)
  {
    x_first <- x@first
    x_last <- invertRleStrand(x@last)

    right_is_first <- which(strand(x_first) == "-")
    idx <- seq_len(length(x))
    idx[right_is_first] <- idx[right_is_first] + length(x)

    ans <- c(x_last, x_first)[idx]
    setNames(ans, names(x))
  }

  left = function(x)
  {
    x_first <- x@first
    x_last <- invertRleStrand(x@last)

    left_is_last <- which(strand(x_first) == "-")
    idx <- seq_len(length(x))
    idx[left_is_last] <- idx[left_is_last] + length(x)

    ans <- c(x_first, x_last)[idx]
    setNames(ans, names(x))
  }

  invertRleStrand = function(x)
  {
    x_strand <- strand(x)
    runValue(x_strand) <- strand(runValue(x_strand) == "+")
    strand(x) <- x_strand
    x
  }

  strand.switch = function(which.strand)
  {
    switch(which.strand,"+" = "-","-" = "+")
  }

}

# preparing rle tracks and coverages
{
  # create transcribed bases rle tracks
  create.transcribed.bases.rle.tracks(bam.files = bam.files,
                                           bam.input.folder = bam.files.path,
                                           rle.out.folder = out.file.path,
                                           prefix = "",
                                           human.chrs = names(chrs.lengths),
                                           human.chrs.lengths = chrs.lengths,
                                           strand.specific = TRUE,
                                           remove.duplicates = FALSE,
                                           size.selection = FALSE,
                                           strand.mode = 2)

  # create cov_i coverages
  coverage.initial.anno.df = data.frame(coverage.initial.anno)
  colnames(coverage.initial.anno.df)[1] = "chr"
  rownames(coverage.initial.anno.df) = coverage.initial.anno.df$ID

  calculate.counts(anno = coverage.initial.anno.df,
                   rle.location = paste0(out.file.path, "/transcribed.bases.rle.tracks/"),
                   bam.files = bam.files,
                   rle.prefix = "",
                   out.folder = out.file.path,
                   file.name = "coverage.initial",
                   antisense.bias.ratio.location = antisense.bias.ratio,
                   chrs = names(chrs.lengths),
                   chrs.lengths = chrs.lengths,
                   antisense = T,
                   find.antisense = F)

  # create cov_f coverages
  coverage.final.anno.df = data.frame(coverage.final.anno)
  colnames(coverage.final.anno.df)[1] = "chr"
  rownames(coverage.final.anno.df) = coverage.final.anno.df$ID

  calculate.counts(anno = coverage.final.anno.df,
                   rle.location = paste0(out.file.path, "/transcribed.bases.rle.tracks/"),
                   bam.files = bam.files,
                   rle.prefix = "",
                   out.folder = out.file.path,
                   file.name = "coverage.final",
                   antisense.bias.ratio.location = antisense.bias.ratio,
                   chrs = names(chrs.lengths),
                   chrs.lengths = chrs.lengths,
                   antisense = T,
                   find.antisense = F)
 }

# drop-off probability estimation
{
  # load coverages
  cov_i = get(load(paste0(out.file.path, "/coverage.initial.antisense.corrected.RData")))
  cov_f = get(load(paste0(out.file.path, "/coverage.final.antisense.corrected.RData")))

  index = rownames(cov_i)
  
  # dividing by respective feature lengths
  cov_i = cov_i[index, ]/coverage.initial.anno.df[index, ]$width
  cov_f = cov_f[index, ]/coverage.final.anno.df[index, ]$width
  gene.lengths = width(anno)

  # ratio of coverages
  rate.ratio = cov_f/cov_i
  dim(rate.ratio)

  # dropoff estimation
  dropoff.probability = 1 - exp((log(rate.ratio))/(gene.lengths/2 - 1))
  dropoff.probability = data.frame(dropoff.probability)
  dropoff.probability[dropoff.probability <= 0 | dropoff.probability >= 1] = NA # setting NA values for genes with cov_i > cov_f
  dropoff.probability = dropoff.probability[-which(is.na(rowSums(dropoff.probability))), ] # subsetting for genes with cov_i > cov_f
  
  # sample averages
  dmso1h = c(1,2) # column numbers, please check and change accordingly
  dmso4h = c(3,4)
  dtag1h = c(5,6)
  dtag4h = c(7,8)
  dropoff.probability["DMSO1h"] = rowMeans(dropoff.probability[, dmso1h])
  dropoff.probability["DMSO4h"] = rowMeans(dropoff.probability[, dmso4h])
  dropoff.probability["DTAG1h"] = rowMeans(dropoff.probability[, dtag1h])
  dropoff.probability["DTAG4h"] = rowMeans(dropoff.probability[, dtag4h])
  
  save(dropoff.probability, file = paste0(out.file.path, "/dropoff.probability.RData"))
  
  # plot (base R plot)
  boxplot(dropoff.probability[, c("DMSO1h", "DMSO4h", "DTAG1h", "DTAG4h")], outline = F, notch = T)
}
  
# figure 3I (using ggplot2)
{
  drop_prob = dropoff.probability[, c(9:12)] %>% tibble %>% 
    gather(., "DMSO1h", "DMSO4h", "DTAG1h", "DTAG4h", key = "sample", value = "drop-off_p") %>%
    mutate(treatment_time = case_when(grepl(sample, pattern = "1h") ~ "1 h",
                                      grepl(sample, pattern = "4h") ~ "4 h",)) %>%
    mutate(Color = case_when(grepl(sample, pattern = "DMSO") ~ "Control",
                             grepl(sample, pattern = "DTAG") ~ paste0(treatment_time," dSSRP1"),))
  drop_prob$Color %<>% factor(., levels = c("Control", "1 h dSSRP1", "4 h dSSRP1"))
  drop_prob$sample %<>% factor(., levels = c("DMSO1h", "DTAG1h", "DMSO4h", "DTAG4h")) 
  FCT.colors = c("blue", "darkblue")
  
  ggplot(data = drop_prob, mapping = aes(x=sample,y=`drop-off_p`,group = sample,fill=Color))+
    geom_boxplot(outlier.alpha = 0.2,outlier.size = 0.3, notch = TRUE, show.legend = TRUE, color = "#000000")+
    stat_compare_means(method = "wilcox.test",comparisons = list(c("DMSO1h", "DTAG1h"), c("DMSO4h", "DTAG4h")), mapping=aes(label=..p.signif..),
                       method.args = list(alternative = "two.sided"),
                       paired = FALSE) +
    scale_fill_manual(values = c("#999999", FCT.colors[1],FCT.colors[2]), name = "")+
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
    labs(title = "TT-seq: drop-off probability [n = 3190]",x="", y = "Probabilty of drop-off")+
    theme_classic(base_family="Helvetica", base_size = 6, base_line_size = 0.5, base_rect_size = 0.5)+
    theme( 
      legend.position = "right",
      legend.direction = "vertical",
      panel.spacing = unit(0.25,units = "cm"),
      plot.title = element_text(hjust = 0.5,color="#000000", size=8),
      strip.text = element_text(color="#000000", size=6),
      strip.background = element_blank(), 
      axis.line = element_line(color="#000000"),
      axis.text = element_text(color="#000000", size=6),
      axis.title = element_text(color="#000000", size=6),
      axis.ticks = element_line(color="#000000"),
      legend.title=element_text(size=6),
      legend.text=element_text(size=6),
      legend.key.size = grid::unit(3,"mm")) 
  
  ggsave(filename = "fig3I.pdf",path = file.path(out.file.path) ,device = "pdf", height=60,width = 60,dpi=300,units="mm",useDingbats=FALSE, family="Helvetica") 
}
