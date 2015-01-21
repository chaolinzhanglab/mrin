suppressMessages(library(getopt, quietly=T))

InputCheck <- function (file)
{
  if (!file.exists (file))
  {
    stop("Input file not exist.") 
  }
	
  dat <- read.table (file, sep="\t", header=T);
  file.mat <- as.matrix (dat[,-1]);
  rownames(file.mat) <- dat[,1];
  #CZ: what if there are duplicates in the gene ID (we do have from the previous step)
  #this needs to be handled properly, here or somewhere else	

  return (file.mat) 

#   if (is.data.frame (file.mat))
#   {
#     return (file.mat)
#   }
#   else
#   {
#     stop ("Wrong Input File Format!")
#   }
}


GetmKS <- function (x)
{
  x.1 <- x[which(!is.na(x))]

  for (i in 1:length(x))
  {
    if(!is.na(x[i]))
    {
      x[i] <- x[i]-median(x.1)  
    }
  }
  x
}

#CZ: the function above is inefficient
GetmKS2 <- function (x)
{
	x = x - median(x, na.rm=T);
	x
}



rpkm.cutoff2 <- function (x)
{
  x.bak <- x 
  x[which(x.bak < 2)] <- 0
  x[which(x.bak >= 2)] <- 1
  x
}

#CZ: whether filtering should be done can be judged from whether gene.matrix is provided

KSFilter  <- function (ks.matrix,  gene.matrix, sd.prop = 0.05)
{
  if ( is.na(gene.matrix))
  {
    mKS <- t(apply(ks.matrix, 1, GetmKS));
    return(mKS);
  }else
  {
    ### Check the order of genes and samples for two matrices ###
    gene.idx <- match(rownames(ks.matrix), rownames(gene.matrix))
    sample.idx <- match(colnames(ks.matrix), colnames(gene.matrix))
    
    gene.matrix.bak <- gene.matrix
    gene.matrix <- gene.matrix[gene.idx, sample.idx]
    
    ### Filter top sd.prop% genes by SD of KS statistics ###
    NA.dis <- apply(ks.matrix, 1, function(x) length(which(is.na(x))))
    #na.num <- 300
	#CZ: this is not good, the line below is just a work around
	na.num <- dim(ks.matrix)[2];

    ### Plot the distribution of NA number of KS statistics and decide the number###
    na.idx <- which(NA.dis < na.num)
    
    ks.matrix <- ks.matrix[na.idx, ]
    gene.matrix<- gene.matrix[na.idx, ]
    
	#CZ: the following line can be simplified by using sd(x, na.rm=T)
    ks.sd <- apply(ks.matrix, 1, function(x) {x <- x[which(!is.na(x))]; return(sd(x))})
    sd.idx <-  order(ks.sd)[floor(length(ks.sd)*sd.prop) : length(ks.sd)]
    
    ### Filter genes that have strong 5' bias ###
	
	#CZ: the following line can be simplified by using median(x, na.rm=T)
    ks.median <- apply(ks.matrix, 1, function(x) {x <- x[which(!is.na(x))]; return(median(x))})
    meidan.cutoff <- 0
    
    include.idx <- intersect(which(ks.median > 0), sd.idx)
    ks.matrix.filtered <- ks.matrix[include.idx, ]
    gene.matrix.filtered <- gene.matrix[include.idx, ]
    
  
    ### Calculate the mKS matrix ###
    mKS <- t(apply(ks.matrix.filtered, 1, GetmKS))
    
    ### Filter mKS by gene expression values with cutoff 2 ###
    rpkm.indicate <- apply(gene.matrix.filtered, 2, rpkm.cutoff2)
    mKS.1 <- mKS * rpkm.indicate
    
    ### Filter mKS by number of mks estimated across samples ### 
    mks.est.num <- apply(mKS.1, 1, function(x) length(which(!is.na(x))))
    include.idx.1 <- which( mks.est.num >= floor(dim(gene.matrix)[2] / 2))
    
    mKS.2 <- mKS.1[include.idx.1,]
  
    return(mKS.2)
  }
}


Cal.mRIN <- function (mKS, ks.matrix)
{
  mRIN <- apply(mKS, 2, function(x) -mean(x[which(!is.na(x))]))
  names(mRIN) <- colnames(ks.matrix)
  mRIN
}

trunc.norm.par <- function(x, alpha)
{
  Z <- 1 - pnorm(alpha)
  
  sigma <- var(x)/(1 + alpha*dnorm(alpha)/Z - dnorm(alpha)^2/Z^2)
  sigma <- sqrt(sigma)
  
  mu <- mean(x) - dnorm(alpha)*sigma/Z
  
  list(mu=mu, sigma=sigma)  

}


KS.stat <- function(par, x, alpha)
{
  #CZ: this is slow and needs to be improved later
  set.seed(84)
  
  est <- rnorm(1000000, par$mu, par$sigma)
  x.prop <- length(x [x > alpha])/length(x)
  y.idx <- order(est, decreasing=T)[1 : floor(length(est) * x.prop)]
  est <- est[y.idx]
  
  ks.stat <- ks.test(x[x > alpha], est)$statistic
  ks.stat
}


Calculate.pvalue <- function(mRIN, verbose=F)
{
  alpha <- seq(-0.05, 0.05, 0.0005)
  
  mu <- alpha
  sigma <- alpha
  dis <- alpha
  
  for (i in 1:length(alpha))
  {

	if (verbose)
	{
		if (T) {cat ('alpha=', alpha[i], '...\n')};
	}
    par <- trunc.norm.par (mRIN [mRIN > alpha[i]],alpha[i])
    mu[i] <- par$mu;
    sigma[i] <- par$sigma;
    dis[i] <- KS.stat (par, mRIN, alpha[i]);
  }
  
  ### Choose best alpha ###
  
  alpha.best <- alpha[which.min(dis)]
  par.best <- trunc.norm.par (mRIN [mRIN > alpha.best],alpha.best)
  
  p.value <- pnorm(mRIN, par.best$mu, par.best$sigma)
  Zscore <- (mRIN - par.best$mu) / par.best$sigma
  
  list(Zscore = Zscore, Pvalue = p.value, alpha.best=alpha.best, mu=par.best$mu, sigma=par.best$sigma)
}

Cal.GIS <- function (mRIN, mKS)
{
  GIS <- NULL
  for (i in 1:dim(mKS)[1])
  {
    idx <- which(!is.na(mKS[i,]))
    GIS <- c(GIS,cor(mKS[i,idx] ,mRIN[idx]))
  }
  names(GIS) <- rownames(mKS)
  GIS  
}


#argument mask: 0 - no argument, 1 - required argument, 2 - optional argument
optionSpec = matrix(c(
    'ksfile',    'k',    1, "character",
    'expfile',   'x',    2, "character",
    'mrinfile',  'm',    1, "character",
    'gisfile',   'G',    1, "character",
    'verbose',   'v',    2, "integer",
    'help'   ,   'h',    0, "logical"
    ), byrow=TRUE, ncol=4);

opt = getopt(optionSpec);


#set some reasonable defaults for the options that are needed,
#but were not specified.

#default parameters
gene.exp.file=NA;
verbose = 0;

ks.file = opt$ksfile;
mrin.file = opt$mrinfile;
gis.file = opt$gisfile;


#set optional arguments
if (!is.null(opt$expfile)) {gene.exp.file = opt$expfile}
if (!is.null(opt$verbose)) {verbose = opt$verbose}

if ( !is.null(opt$help) |is.null(opt$ksfile) | is.null(opt$mrinfile)| is.null(opt$gisfile))
{
    #cat(getopt(optionSpec, usage=TRUE));
    cat (
        'calculate mRIN and GIS from KS matrix\n',
        'Usage: Rscript ', get_Rscript_filename(),"\n",
        '[required]\n',
        ' -k, --ksfile     [string]: input file with KS matrix\n',
        ' -m, --mrinfile   [string]: output file for mRIN\n',
        ' -G, --gisfile    [string]: output file for GIS\n',
        '[options]\n',
        ' -x, --expfile    [string]: input file with gene RPKM values\n',
        ' -v, --verbose            : verbose mode\n',
        ' -h, --help               : print usage\n'
    );
    q(status=1);
}


### Check Input file ###
if (verbose) {cat ('read KS matrix file ...\n');}
ks <- InputCheck (ks.file)

gene.exp <- NA;
if (!is.na(gene.exp.file))
{
	if (verbose) {cat ('read gene expresion matrix file ...\n');}
	gene.exp <- InputCheck (gene.exp.file)
}

### Calculate mRIN for each sample ###
if (verbose) {cat ('filter KS matrix ...\n');}
mKS <- KSFilter (ks, gene.exp)

if (verbose) {cat ('calculate sample mRINs ...\n');}
mRIN <- Cal.mRIN (mKS,ks)

if (verbose) {cat ('estimate p-values ...\n');}
mRIN.stat <- Calculate.pvalue (mRIN, verbose)
Pvalue <- mRIN.stat$Pvalue
Zscore <- mRIN.stat$Zscore

if (verbose) {cat ('calculate GIS ...\n');}
GIS <- Cal.GIS (mRIN, mKS)

if (verbose) {cat ('save output ...\n');}
write.table(cbind(mRIN,Zscore, Pvalue), mrin.file, sep="\t", quote=FALSE, col.names=T, row.names=T)
write.table(cbind(names(GIS), GIS), gis.file, col.names=c("gene.symbol","GIS"), quote=F, sep="\t", row.names=F)


