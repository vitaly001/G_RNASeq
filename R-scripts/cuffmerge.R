setwd("/Volumes/HD2/G_RNASeq")
samples = read.csv("samples.csv", stringsAsFactors=FALSE)

setwd("/Volumes/HD2/G_RNASeq/Cufflinks")

gf = "/Volumes/HD2/UCSC/mm10/Annotation/Archives/archive-2015-07-17-14-33-26/Genes/genes.gtf"
for(i in seq_len(nrow(samples))) {
        lib = samples$conditions[i]
#        print(lib)
        cuffFile = file.path(lib, "transcripts.gtf")
        cat(cuffFile, file = "assem_GTF.txt", collapse = '\n', append = TRUE)
}


system(paste0("nohup cuffmerge -g ", gf ," -p 6 ", " -o", " /Volumes/HD2/G_RNASeq/Cuffmerge"," /Volumes/HD2/G_RNASeq/Cufflinks/assem_GTF.txt", " >& com.cuffmerge.log &"))

