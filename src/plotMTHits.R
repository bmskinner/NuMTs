# Read mt sam data
library(dplyr)
library(ggplot2)
# source("https://bioconductor.org/biocLite.R")
# biocLite("GenomicRanges")
# biocLite("ggbio")
library(GenomicRanges)
library(ggbio)

files = c("SRR7142641.filt.MT.sam", "SRR7142642.filt.MT.sam", "SRR7142643.filt.MT.sam" )

# Pig genome data
chrnames   = c(1:18, "X", "Y")
chrlengths = c(274330532, 151935994, 132848913, 130910915, 104526007, 170843587, 121844099, 138966237, 139512083, 69359453, 79169978, 61602749, 208334590, 141755446, 140412725, 79944280, 63494081, 55982971, 125939595, 43547828)
names(chrlengths) = chrnames
refgenome = GRanges(chrnames, IRanges(0, chrlengths)))

readSamFile = function(file){
  data = read.csv(file, skip=1, header=F, stringsAsFactors = F, sep="\t")
  locs = data %>% select(V3, V4, V7, V8, V10) %>% filter(V3=="MT"| V7=="MT")
  
  r1 = locs %>% filter(V3=="MT") %>% select("MT_pos" = V4,"Chr"=V7,"Chr_pos"=V8, "Seq"=V10  )
  r2 = locs %>% filter(V3!="MT") %>% select("MT_pos" = V8,"Chr"=V3,"Chr_pos"=V4, "Seq"=V10  )
  
  mge = rbind(r1, r2) %>% mutate(TLen = nchar(Seq)) %>% select(-Seq)
  mge$File = file
  
  seqs = reduce(with(mge, GRanges(Chr, IRanges(Chr_pos, Chr_pos+TLen))))
  seqlevels(seqs) = chrnames
  seqlengths(seqs) = chrlengths
  
  seqs
}

seqs = lapply(files, readSamFile)

sub2 = Reduce(subsetByOverlaps, seqs)
seqlevels(sub2) = chrnames
seqlengths(sub2) = chrlengths

ggbio() + circle(sub2, geom = "ideo", fill = "gray70") +
  circle(sub2, geom = "rect", color = "red", size=2) +
  circle(sub2, geom = "scale", size = 2)+
  circle(sub2, geom = "text", aes(label = seqnames), vjust = 0, size = 3)
