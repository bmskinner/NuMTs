# This script parses the pig SRA dataset to download the SRA data
# for animals with gDNA and defined breeds

# cat("Loading R packages ...\n")
suppressMessages(library(dplyr))
suppressMessages(library(stringr))

base_dir     = commandArgs(TRUE)[1]
srrAccession = commandArgs(TRUE)[2]
breed        = commandArgs(TRUE)[3]

# #' This function loads the external functions file.
# #' @title Load external functions
# #' @export
# loadFunctionsFile = function(){
#   cat("Loading functions ...\n")
#   file.arg.name = "--file="
#   script.name   = sub(file.arg.name, "", commandArgs()[grep(file.arg.name, commandArgs())])
#   script.folder = dirname(script.name)
#   script.to.load = paste(sep="/", script.folder, "functions.r")
#   source(script.to.load)
# }

# loadFunctionsFile()

# base_dir = pathExistsOrQuit(base_dir, "Base analysis directory")


# Read the accession table
# samples = read.csv(paste0(base_dir,"Filtered_samples.csv"), header=T, stringsAsFactors = F, sep="\t")

# Filter the SRA data to remove ambiguous breeds and samples without enough data
# filterSraRuns = function(data){
#   ignoreBreeds = c("", "not collected", "NS", "X bred", "pigs", "min", "F1", "missing","Feral pig", "Min pig", "Cross-bred", "northeast wild", "domestic", "Chinese domestic pigs", "Jiangsu pig breeds", "White Composite", "Asian domestic pig", "Asia wild" )
  
#   filt = data %>% 
#     filter(Assay_Type=="WGS", 
#                          !(breed %in% ignoreBreeds), 
#                          LibraryLayout=="PAIRED",
#                          MBases>100) %>% 
#     mutate(breed = gsub("Durocs", "Duroc", breed)) %>%
#     mutate(breed = gsub("LargeWhite", "Large White", breed))
  
#   summ = filt %>% group_by(breed) %>% mutate(nSamples = n()) %>% filter(nSamples>3) %>% select(breed, nSamples) %>% distinct()

#   samples = filt %>% group_by(breed) %>% mutate(nSamples = n()) %>% 
#     select(breed, BioProject, BioSample, Experiment, MBases, SRA_Sample, Run, sex) %>%
#     arrange(breed)
#   samples
# }

# samples = filterSraRuns(data)

# write.table(samples, file=paste0(base_dir, "Filtered_samples.csv"), col.names=T, row.names=F, sep="\t")

# Download FASTQ files directly from the EBI
downloadSrrFromEbi = function(srrAccession, breed){

  baseDir = paste0(base_dir, "fastq/", breed, "/")
  
  # Make a breed directory if needed
  dir.create(baseDir)

  dir1 = stringr::str_extract(srrAccession, "^(\\w{3}\\d{3})")
  endDigit = stringr::str_extract(srrAccession, "(\\d)$")
  isDir2 = str_detect(srrAccession, "^\\w{3}\\d{7}$")
  
   # only used if >6 digits in accession
  dir2 = ifelse(isDir2, paste0("00", endDigit, "/"), "")
    
  fq1 = "_1.fastq.gz"
  fq2 = "_2.fastq.gz"
  
  urlBase = paste0("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/", dir1, "/", dir2, srrAccession, "/", srrAccession)
  
  url1 = paste0(urlBase, fq1)
  url2 = paste0(urlBase, fq2)
  
  out1 = paste0(baseDir, srrAccession, fq1)
  out2 = paste0(baseDir, srrAccession, fq2)

  download.fastq = function(url, out){
    if(!file.exists(out)){
        cat(breed, ": ", srrAccession, ": Downloading: ", url, "\n")
        tryCatch(
            download.file(url, method = "libcurl", destfile=out, quiet=T )
            ,error=function(e){
              cat("Unable to download ", srrAccession, ":", e,"\n")
            }
        )
      } else{
        cat("Skipping existing:", out,"\n")
      }
  }

  download.fastq(url1, out1)
  download.fastq(url2, out2)  
}

# download all desired runs
downloadSrrFromEbi(srrAccession, breed)