options(warn = -1)
#if(!require('Matrix')){install.packages('Matrix')}
args <- commandArgs(trailingOnly = TRUE)
inputfile <- as.character(args[1])
outfile <- as.character(args[2])
t1<-Sys.time()
rawcount <- read.table(inputfile, header=T)
#rawcountbak <- rawcount
colnames(rawcount) <- paste0(colnames(rawcount),"-1")
rawcount <- Matrix::Matrix(as.matrix( rawcount ), sparse = TRUE)
t2<-Sys.time()
rawcount_fix <- rawcount
allnames <- rownames(rawcount)
overlapid <- grep(",",allnames)
for( id in overlapid ){
	idn <- unlist(strsplit(allnames[id],split=","))
	rawcount_fix[idn,] <- apply(rawcount_fix[c(idn,allnames[id]),],2,function(x){ idl <- length(x); y <- (x[-idl]+x[idl]/(idl-1)) ; return(y) })
}
rawcount_fix <- rawcount_fix[-overlapid,]
rawcount_fix <- round(rawcount_fix,0.01)
write.table( as.matrix(rawcount_fix), file = outfile, quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")
system(paste0("gzip -f ",outfile))
