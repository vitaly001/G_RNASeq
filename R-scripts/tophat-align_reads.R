#Align the reads (using tophat2) to reference genome
#Using R string manipulation, construct the Unix commands to call tophat2
getwd()
setwd("/Volumes/HD2/G_RNASeq/")
samples = read.csv("sampletab.csv", stringsAsFactors=FALSE)


gf = "/Volumes/HD2/UCSC/mm10/Annotation/Archives/archive-2015-07-17-14-33-26/Genes/genes.gtf"
bowind = "/Volumes/HD2/UCSC/mm10/Sequence/Bowtie2Index/genome"
cmd = with(samples,
           paste("tophat -G", gf, "-p 5 -o", conditions, bowind,
                 fastq1))

align = do.call(paste, c(as.list(cmd), sep =' && ')) #concatenate all commands in vector separated with shell control operator

align # print the consecutive commands && : Used to build AND lists, it allows you to run one command only if another exited successfully. 
system(align)  # invoke commands
# system(cmd) # invoke command
#for (i in cmd) {system(i)}  # invoke commands using loop (different method)
