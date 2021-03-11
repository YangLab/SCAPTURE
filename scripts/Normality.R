options(warn = -1)
#if(!require('nortest')){install.packages('nortest')}
#if(!require('diptest')){install.packages('diptest')}
args <- commandArgs(trailingOnly = TRUE)
genepeak <- as.character(args[1])
genecov <- as.character(args[2])
peakwidth <- as.numeric(args[3])/2
cov <- read.table(genecov, header=F, stringsAsFactors = F)
peak <- read.table(genepeak)

pvalue.nomal = c()
pvalue.symmetry = c()
pvalue.bimodal = c()

for(j in 1:nrow(peak)){
    s <- as.numeric(peak[j,3])
    e <- as.numeric(peak[j,4])
    n <- which.max(cov$V5[s:e]) + s
    d <- min(c(n-s, e-n))
    s <- n-d
    e <- n+d
    if( d <= peakwidth/10 ){s = as.numeric(peak[j,3]);e = as.numeric(peak[j,4]);d = (e-s)/2;}
    
    #decrease number of data point
    magn10 <- 1
    magn <- log(mean(cov$V5[s:e],na.rm = T),10)
    if(magn >= 3){ magn10 <- 10^(floor(magn) - 3) }

    dat <- c()
    dat_l <- c()
    dat_r <- c()
    dat_l <- round(cov$V5[s:n]/magn10,0)
    dat_r <- round(cov$V5[n:e]/magn10,0)
    for(k in s:e){ r <- round(cov$V5[k]/magn10,0); r <- ifelse(is.na(r),0,r); dat = c(dat, rep(k-s-d, r) ) }

    ifnomal <- c() #recommend 1e-4 as cutoff
    ifreplace <- FALSE
    if(length(dat) <= 200){ ifreplace = TRUE }
    for(pn in 1:100){ 
        p <- try(nortest::ad.test(dat[sample(1:length(dat), 200, replace = ifreplace)])$p.value, silent=TRUE)
        if(!is.numeric(p)){p <- 0}
        ifnomal = c(ifnomal, p) }
    
    #symmetry distribution test
    ifsymmetry <- c() #recommend 1e-4 as cutoff
    ifreplace <- FALSE
    if(length(dat_l) <= 15 | length(dat_r) <= 15){ ifreplace = TRUE }
    for(pn in 1:100){
        p <- try(t.test(dat_l[sample(1:length(dat_l), 15, replace = ifreplace)],dat_r[sample(1:length(dat_r),15,replace = FALSE)])$p.value, silent=TRUE)
        if(!is.numeric(p)){p <- 0}
        ifsymmetry = c(ifsymmetry, p) }
    
    #bimodal distribution test
    ifbimodal <- c()
    ifreplace <- FALSE
    if(length(dat) <= 5000){ ifreplace = TRUE }
    ifbimodal <-  try(diptest::dip.test(dat[sample(1:length(dat),5000)])$p.value, silent=TRUE)
    if(!is.numeric(ifbimodal)){ifbimodal <- 1}  

    #output pvalue
    pvalue.nomal <- c(pvalue.nomal, median(ifnomal))
    pvalue.symmetry <- c(pvalue.symmetry, median(ifsymmetry))
    pvalue.bimodal <- c(pvalue.bimodal, ifbimodal)
}
write.table(cbind(peak,pvalue.nomal,pvalue.symmetry,pvalue.bimodal),file = sub(".txt",".Pvalue.txt",genepeak), col.names = F, row.names = F, quote = F, sep = "\t")

