# Plot MT reads using bam file

library(dplyr)
library(ggplot2)
# source("https://bioconductor.org/biocLite.R")
# biocLite("GenomicRanges")
# biocLite("ggbio")
library(GenomicRanges)
library(GenomicAlignments)
library(ggbio)

# Pig genome data
makeRefGenome = function(){
  chrnames   = c(1:18, "X", "Y", "MT")
  chrlengths = c(274330532, 151935994, 132848913, 130910915, 104526007, 170843587, 121844099, 138966237, 139512083, 69359453, 79169978, 61602749, 208334590, 141755446, 140412725, 79944280, 63494081, 55982971, 125939595, 43547828, 16613)
  names(chrlengths) = chrnames
  gr = GRanges(chrnames, IRanges(1, chrlengths), seqlengths = chrlengths)
  gr = sortSeqlevels(gr)
  gr = sort(gr)
  gr
}

genome.bam = makeRefGenome()

#' Get the reads from a bam file
#'
#' @param file the file to read
#'
#' @return a GRanges object containing the reads
#' @export
getReadsFromBam = function(file){
  # cat("Opening:", file,"\n")
  
  # Read the bam file
  result = tryCatch({
    bamFile = BamFile(file, asMates=TRUE, index = paste0(file, ".bai"))
    gapairs = readGAlignmentsList( bamFile)
    
    # Calculate coverage and cut at <min.coverage>
    covg = coverage(gapairs)
    min.coverage = mean(mean(covg))
    peaks = XVector::slice(covg, lower=min.coverage)
    ranges = as(peaks, "GRanges") 
    red = reduce(ranges)
    gr = sortSeqlevels(red)
    gr = sort(gr)
    gr
    
  }, error=function(error){
    message(error)
    cat("\n")
    NA
  })
  result
}


#' Reduce all Granges in a list 
#'
#' @param list 
#'
#' @return
#' @export
#'
#' @examples
uniqueGranges = function(list){
  result = list[[1]]
  
  for(ranges in list){
    
    result = c(result, ranges)
  }
  result = reduce(result)
  cat("Found", length(result),"unique ranges\n")
  result
}



findSharedInsertions = function(breed){
  
  cat("Running on", breed,"\n")
  
  ## Load data files and reference genome
  dir = paste0("Y:/bms41/Pigs/NuMTs/filtered_reads/",breed,"/")
  files = list.files(dir, pattern=".*.MT.bam$", all.files=F, full.names=T)
  
  # Create a list of all ranges in the bam files
  all.seqs = lapply(files, getReadsFromBam)
  all.seqs = all.seqs[!is.na(all.seqs)]
  cat(length(all.seqs),"files opened\n")

  cat("Getting unique ranges\n")
  unique.ranges = uniqueGranges(all.seqs)
  
  if(length(unique.ranges)==0){
    return(NA)
  }
  
  presentInAtLeast = function(list, unique.ranges, percent){
    
    min.overlaps = (length(list) * percent)/100
    cat("Looking for at least", min.overlaps, "overlaps for each range\n")
    
    countOverlappingGRanges = function(test.range){countOverlaps(unique.ranges, test.range ) }
    
    # Find the overlapping ranges
    overlaps = lapply(list, countOverlappingGRanges)
    
    processList = function(list){
      count = rep(0, length(list[[1]]))
      
      for(l in list){
        count = mapply(sum, count, l)
      }
      count
    }
    
    # count the overlapping ranges
    totals = processList(overlaps)
    
    # filter by count
    filterCounts = function(count){ count>min.overlaps }
    
    keeps = filterCounts(totals)
    
    # Create a result variable by subsetting the first range in the list
    result = unique.ranges[1]
    
    for(i in 1:length(unique.ranges)){
      if(keeps[[i]]){
        result = c(result, unique.ranges[i])
      }
    }
    # First row was a placeholder, remove unless in the keep list
    if(!keeps[[1]]){ result = result[2:length(result)]  }
    cat("Found", length(result),"ranges in at least", percent,"% of the samples\n")
    reduce(result)
  }
  
  shared.ranges = presentInAtLeast(all.seqs, unique.ranges, 60)
  seqlengths(shared.ranges) = seqlengths(genome.bam)
  shared.ranges = sortSeqlevels(shared.ranges)
  shared.ranges = sort(shared.ranges)
  shared.ranges$breed = breed
  shared.ranges
}




plotInserts = function(range.list){
  g = ggbio() + circle(genome.bam, geom = "ideo", fill = "gray70")
  for(s in range.list){
    # s = applyChrSizes(s)
    cat("Drawing", unique(s$breed), "\n")
    g = g + circle(s, geom = "rect", color = "steelblue", size=1)
    
  }
  g = g +
    # circle(shared.ranges, geom = "rect", color = "red", size=2) +
    circle(genome.bam, geom = "scale", size = 2)+
    circle(genome.bam, geom = "text", aes(label = seqnames), vjust = 0, size = 3)
  g
  ggsave("Circle.png", g, dpi=300)
}


writeInsertTables = function(range.list){
  for(s in range.list){
    cat("Writing", unique(s$breed), "\n")
    df = data.frame(chr=seqnames(s),
                           start=start(s),
                           end=end(s))
    
    file.name = paste0("Hits.", unique(s$breed), ".txt")
    write.csv(df, file=file.name, row.names = F, col.names = T)
  }
    
}

breeds = c("Duroc", "Berkshire", "Landrace")
inserts = lapply(breeds, findSharedInsertions)
inserts = inserts[!is.na(inserts)]
plotInserts(inserts)
writeInsertTables(inserts)
