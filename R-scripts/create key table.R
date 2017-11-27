#set working directory
setwd("/Volumes/HD2/G_RNAseq")

#create an character vector with path to all csv files
LF = list.files(recursive=T, pattern = paste('*.csv', sep=''))

#constructs and executes a function call to read all csv files in all subfolders and rbind the rows together

dat = do.call(rbind, lapply(LF, read.csv))

#write resulted dataframe as csv file

write.csv(dat, file="sampletab.csv", row.names=F)
