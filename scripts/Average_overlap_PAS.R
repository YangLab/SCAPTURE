options(warn = -1)
#if(!require('Matrix')){install.packages('Matrix')}
#if(!require('parallel')){install.packages('parallel')}
args <- commandArgs(trailingOnly = TRUE)
inputfile <- as.character(args[1])
outfile <- as.character(args[2])
rawcount <- read.table(inputfile, header=T)
colnames(rawcount) <- paste0(colnames(rawcount),"-1")
rawcount <- Matrix::Matrix(as.matrix( rawcount ), sparse = TRUE)
rawcount_fix <- rawcount
allnames <- rownames(rawcount)
overlapid <- grep(",",allnames)
if(length(overlapid) == 0){
	system(paste0("cp ",inputfile," ",outfile,".gz"))
	quit(save="no")
}
overlapname <- allnames[overlapid]
multiplethreads <- TRUE

if(multiplethreads == FALSE){
for( id in overlapid ){
	idn <- unlist(strsplit(allnames[id],split=","))
	rawcount_fix[idn,] <- apply(rawcount_fix[c(idn,allnames[id]),],2,function(x){ idl <- length(x); y <- (x[-idl]+x[idl]/(idl-1)) ; return(y) })
}
rawcount_fix <- rawcount_fix[-overlapid,]
rawcount_fix <- round(rawcount_fix,0.01)
}

if(multiplethreads == TRUE){
mc <- getOption("mc.cores", 4)
PAScollaspe <- parallel::mclapply(overlapid,
    function(id){
        idn <- unlist(strsplit(allnames[id],split=","))
        fix <- apply(rawcount_fix[c(idn,allnames[id]),],2,function(x){ idl <- length(x); y <- (x[-idl]+x[idl]/(idl-1)) ; return(y) })
        return(fix)},
    mc.cores = mc)

PAScollaspeMatrix <- c()
for( i in 1:length(PAScollaspe)){
    PAScollaspeMatrix <- rbind(PAScollaspeMatrix,PAScollaspe[[i]])
}
PAScollaspeMatrix <- Matrix::Matrix(as.matrix( PAScollaspeMatrix ), sparse = TRUE)

rmid <- match(c(rownames(PAScollaspeMatrix),overlapname),rownames(rawcount_fix))
rawcount_fix <- rawcount_fix[-rmid,]
rawcount_fix <- rbind(rawcount_fix,PAScollaspeMatrix)
rawcount_fix <- round(rawcount_fix,0.01)
}

write.table( as.matrix(rawcount_fix), file = outfile, quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")
system(paste0("gzip -f ",outfile))
