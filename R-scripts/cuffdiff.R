setwd("/Volumes/HD3/NGS/G_RNASeq")
samples = read.csv("samples.csv", stringsAsFactors=FALSE)
genome = "/Volumes/HD3/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa"
gf = "/Volumes/HD3/NGS/G_RNASeq/Cuffmerge/merged.gtf"
rm(testfiles)
testfiles = ''
cuff_labels = ''

for(i in seq_len(nrow(samples))) {
        lib = samples$conditions[i]
        bamFile = file.path(lib, "accepted_hits.bam")
        cuff_labels = paste0(cuff_labels, lib, sep = ",")
        testfiles = paste0(testfiles, bamFile, sep = ' ')
        
        
}
# cuffdiff should run with genome.fa file and option -b
print(paste0('cuffdiff -o Cuffdiff -p 6 ', ' -L ', cuff_labels, " -b ", genome, ' -u ', gf, ' ', testfiles, ' >& com.cuffdiff.log &'))

print(testfiles)


#samples in each condition should be separated by comma, between condition should be space only
#system(paste0("nohup cuffmerge -g ", gf ," -p 6 ", " -o", " /Volumes/HD2/ngsmus/Cuffmerge"," /Volumes/HD2/ngsmus/Cufflinks/assem_GTF.txt", " >& com.cuffmerge.log &"))
        
