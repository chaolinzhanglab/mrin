suppressMessages(library(getopt, quietly=T))

InputCheck <- function (file)
{
  if (!file.exists (file))
  {
    stop("Input file not exist.") 
  }
	
  #CZ: the original code has rownames, which is not compatible with output of the previous step
  #and the idea of having a column without row header is weird anyway ...

  dat <- read.table (file, sep="\t", header=T);
  file.mat <- as.matrix (dat[,-1]);
  rownames(file.mat) <- dat[,1];
  #CZ: what if there are duplicates in the gene ID (we do have from the previous step)
  #this needs to be handled properly, here or somewhere else	

  return (file.mat) 
}



GetmKS <- function (x)
{
	x = x - median(x, na.rm=T);
	x
}



rpkm.filter <- function (x, cutoff)
{
  x.bak <- x 
  x[which(x.bak < cutoff)] <- 0
  x[which(x.bak >= cutoff)] <- 1
  x
}


KSFilter  <- function (ks.matrix,  gene.matrix, propSD, fivebias, minRPKM, propEmks)
{
  ### Remove duplicate genes ###
  ks.matrix <- ks.matrix[!duplicated(rownames(ks.matrix)), ]
  if ( is.null(dim(gene.matrix)))
  {
    mKS <- t(apply(ks.matrix, 1, GetmKS));
    return(mKS);
  }else
  {
    ### Check if the RPKM cutoff is set by users ###
    if (is.na(minRPKM))
    {
	stop('RPKM cutoff must be set\n')
    }
    ### Check the order of genes and samples for two matrices ###
    gene.matrix <- gene.matrix [!duplicated(rownames(gene.matrix)), ]
    gene.idx <- match(rownames(ks.matrix), rownames(gene.matrix))
    sample.idx <- match(colnames(ks.matrix), colnames(gene.matrix))
    
    gene.matrix.bak <- gene.matrix
    gene.matrix <- gene.matrix[gene.idx, sample.idx]
    
    ### Filter top sd.prop% genes by SD of KS statistics ###
    if( !is.na(propSD))
    {
      if(propSD > 0 & propSD < 1)
      {
        ks.sd <- apply(ks.matrix, 1, function(x)  sd(x, na.rm=T))
        sd.idx <-  order(ks.sd)[floor(length(ks.sd)*propSD) : length(ks.sd)]
      }else
      {
        stop('Wrong propSD!\n')
      }
      
    }
    
    
    ### Filter genes that have strong 5' bias ###
    if(! is.na(fivebias))
    {
	    ks.median <- apply(ks.matrix, 1, function(x) median(x, na.rm=T))
      ks.median.idx <- which(ks.median > 0)
	  }
    
    
    if( !is.na(propSD) & is.na(fivebias))
    {
      ks.matrix.filtered <- ks.matrix[sd.idx, ]
      gene.matrix.filtered <- gene.matrix[sd.idx, ]
    }else if( is.na(propSD) & !is.na(fivebias))
    {
      ks.matrix.filtered <- ks.matrix[ks.median.idx, ]
      gene.matrix.filtered <- gene.matrix[ks.median.idx, ]
    }else if(!is.na(propSD) & !is.na(fivebias))
    { 
      ks.matrix.filtered <- ks.matrix[intersect(sd.idx, ks.median.idx), ]
      gene.matrix.filtered <- gene.matrix[intersect(sd.idx, ks.median.idx), ]
    }else
    {
      ks.matrix.filtered <- ks.matrix
      gene.matrix.filtered <- gene.matrix
    }
	 
	    
  
    ### Calculate the mKS matrix ###
    mKS <- t(apply(ks.matrix.filtered, 1, GetmKS))
    
    ### Filter mKS by gene expression values with cutoff 2 ###
    if (!is.na(minRPKM))
    {
      rpkm.indicate <- apply(gene.matrix.filtered, 2, function(x) rpkm.filter(x, minRPKM))
      mKS <- mKS * rpkm.indicate
    }
    
    
    ### Filter mKS by number of mks estimated across samples ###
    
    if (!is.na(propEmks))
    {
      mks.est.num <- apply(mKS, 1, function(x) length(which(!is.na(x))))
      include.idx.1 <- which( mks.est.num >= floor(dim(gene.matrix)[2] / 2))
      mKS <- mKS[include.idx.1,]
    }
    
    return(mKS)
  }
}


Cal.mRIN <- function (mKS)
{
  mRIN <- apply(mKS, 2, function(x) -mean(x, na.rm=TRUE))
  #names(mRIN) <- colnames(ks.matrix)
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
  alpha <- seq(-0.05, max(mRIN), 0.0005)
  alpha.idx <- which(sapply(alpha,function(x) length(mRIN[mRIN>x]))>1)
  alpha<- alpha[alpha.idx]
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
    'minRPKM',   'r',    2, "double",
    'propSD',    's',    2, "double",
    'propEmks',  'e',    2, "double",
    'fivebias',  'b',    2, "character",
    'verbose',   'v',    2, "integer",
    'help'   ,   'h',    0, "logical"
    ), byrow=TRUE, ncol=4);

opt = getopt(optionSpec);


#set some reasonable defaults for the options that are needed,
#but were not specified.

#default parameters
gene.exp.file=NA;
verbose = 0;
minRPKM = NA;
propSD = NA;
propEmks = NA;
fivebias = NA;

ks.file = opt$ksfile;
mrin.file = opt$mrinfile;
gis.file = opt$gisfile;


#set optional arguments
if (!is.null(opt$expfile)) {gene.exp.file = opt$expfile}
if (!is.null(opt$verbose)) {verbose = opt$verbose}
if (!is.null(opt$minRPKM)) {minRPKM = opt$minRPKM}
if (!is.null(opt$propSD)) {propSD = opt$propSD}
if (!is.null(opt$propEmks)) {propEmks = opt$propEmks}
if (!is.null(opt$fivebias)) {fivebias = opt$fivebias}

if ( !is.null(opt$help) |is.null(opt$ksfile) | is.null(opt$mrinfile)| is.null(opt$gisfile))
{
    #cat(getopt(optionSpec, usage=TRUE));
    cat (
        'calculate mRIN and GIS from KS matrix\n',
        'Usage: Rscript ', get_Rscript_filename(),"\n",
        '[required]\n',
        ' -k, --ksfile     Input file with KS matrix\n',
        ' -m, --mrinfile   Output file for mRIN\n',
        ' -G, --gisfile    Output file for GIS\n',
        '[options]\n',
        ' -x, --expfile    Input file with gene RPKM values\n',
        ' -r, --minRPKM    Gene RPKM value cutoff for ks statistics filtering\n',
        ' -s, --propSD     Proportion of genes with smallest standard deviation of ks statistic for filtering, number between 0 and 1\n',
        ' -e, --propEmks   Proportion of samples with mKS estimated for ks statistic filtering,number between 0 and 1\n',
        ' -b, --fivebias   Remove genes with strong 5 bias\n',
        ' -v, --verbose    Verbose mode\n',
        ' -h, --help       Print usage\n'
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
if (verbose) {cat ('calculate mKS matrix ...\n');}
mKS <- KSFilter (ks, gene.exp, propSD, fivebias, minRPKM, propEmks)

if (verbose) {cat ('calculate sample mRINs ...\n');}
mRIN <- Cal.mRIN (mKS)

if (verbose) {cat ('estimate p-values ...\n');}
mRIN.stat <- Calculate.pvalue (mRIN, verbose)
Pvalue <- mRIN.stat$Pvalue
Zscore <- mRIN.stat$Zscore

if (verbose) {cat ('calculate GIS ...\n');}
GIS <- Cal.GIS (mRIN, mKS)

if (verbose) {cat ('save output ...\n');}
output.mRIN <- cbind(colnames(ks), mRIN, Zscore, Pvalue)
colnames(output.mRIN) <- c("Sample", "mRIN", "Zcore", "Pvalue")
write.table(output.mRIN, mrin.file, sep="\t", quote=FALSE, col.names=T, row.names=F)
write.table(cbind(names(GIS), GIS), gis.file, col.names=c("Gene symbol","GIS"), quote=F, sep="\t", row.names=F)


