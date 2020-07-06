


#' Learns the parameters of a 2 states HHM
#'
#' @param bw_file Path to the input bigWig file
#' @param seqlen The genomic sizes (in bp) of the training chromosomes
#' @param winsize The bin size to use (in bp). Default: 5000.
#' @param stepsize The step size to use. Each time the window is slided by this value. Default: 5000.
#' @param plot.model Whether to plot the learned model or not. Default: FALSE.
#' @param blacklist A GRanges object containing the coordinates of the blacklist regions to exclude. Default: NULL
#' @param smooth Do you want to smooth the signal first then estimate paramters? Generally useful when the signal is sparse. Default: FALSE.
#' @param w The number of bins used to smooth the signal. Default: 3. Ignored if smooth = FALSE.
#'
#' @return
#'  A \LinkS4Class{hmmspec} object containing the paramters of the learned HHM model.
#' @export
#'
#' @examples
#' 
#' bw_file = "./bigWig_examples/H3K27me3_4C_mm10.sorted.Q10.dedup.sorted.bw"
#' hmmmodel = TrainHMMmodel(bw_file = bw_file,
#'                          seqlen = seqlengths(BSgenome.Mmusculus.UCSC.mm10)[1:4],
#'                          blacklist = mm10_blacklist_gr,
#'                          plot.model = FALSE,
#'                          winsize = winsize,
#'                          step=stepsize,
#'                          smooth = FALSE
#'                          )

TrainHMMmodel <- function(bw_file, 
                          seqlen, 
                          winsize=5000, 
                          stepsize=5000, 
                          plot.model = FALSE,
                          blacklist=NULL, 
                          smooth=FALSE, 
                          w=3){
  
  cat(green("Estimating the initial distribution P0\n"))
  
  res <- estimateInitialDomainsDistribution(bw_file, seqlen,winsize = winsize, stepsize = stepsize,blacklist, smooth, w)
  initial <-  res$P0
  
  P <- matrix(c(0.90,0.1,0.1,0.90),2,2)
  b <- list(mu = res$mu, sigma= res$sigma)
  
  ## Create initial model
  cat(green("defining model\n"))
  f = function (x, j, model)
  {
    ret = dnorm(x, model$parms.emission$mu[j], sqrt(model$parms.emission$sigma[j]))
    ret[is.na(ret)] = 1
    #ret[is.infinite(ret)] = 1
    ret
  }
  
  model <- hmmspec(init = initial, trans = P, parms.emission = b, dens.emission = f)
  
  ##
  
  scores <- c()
  written = FALSE
  for(ch in names(seqlen)){
    region <- GRanges(ch, IRanges(1,seqlen[ch]))
    bins = slidingWindows(region,width = winsize,step = stepsize)[[1]]
    
    if(!is.null(blacklist)){
      if(! written == FALSE) {
        cat(yellow("removing blacklist regions\n"))
        written = TRUE
      }
      bins = subsetByOverlaps(bins,blacklist,invert = T)
    }
    
    scores <- c( scores, calcPeaksSignal(bins, bw_file)$meanScore)
  }
  
  q99 = quantile(scores,0.99)
  if(max(scores)/q99>4){
    scores[scores>q99] = q99
  }
  
  
  if(smooth){
    cat(paste0("Smoothing ",ch," scores"))
    scores = smoothSignal(scores,window = w)
  }
  
  cat(green("Fitting model\n"))
  train <- list(x=scores, N= length(scores))
  bw_fileModel <- hmmfit(train, model, mstep=mstep.norm,maxit = 100)$model
  
  
  if(plot.model){
    
    hist(scores, probability=TRUE, breaks=30, xlab="signal enrichment", ylab="Probability", ylim=c(0,1), main="")
    x = seq(min(scores),max(scores),0.01)
    lines(x, dnorm(x, mean=bw_fileModel$parms.emission$mu[1], sd=sqrt(bw_fileModel$parms.emission$sigma[1])),col="red",lwd=2);
    lines(x, dnorm(x, mean=bw_fileModel$parms.emission$mu[2], sd=sqrt(bw_fileModel$parms.emission$sigma[2])),col="green",lwd=2);
  }
  return(bw_fileModel)
}



#' Estimates the initial paramters of the model
#'
#' @param bw_file Path to the input bigWig file
#' @param seqlen The genomic sizes (in bp) of the training chromosomes
#' @param winsize The bin size to use (in bp). Default: 5000.
#' @param stepsize The step size to use. Each time the window is slided by this value. Default: 5000.
#' @param blacklist A GRanges object containing the coordinates of the blacklist regions to exclude. Default: NULL
#' @param smooth Do you want to smooth the signal first then estimate paramters? Generally useful when the signal is sparse. Default: FALSE.
#' @param w The number of bins used to smooth the signal. Default: 3. Ignored if smooth = FALSE. 
#'
#' @return
#' 
#' Resturns a list containing the initial paramters of the model.
#' 
#' \itemize{
#'  \item \code{P0} : The initial probability of observing each state (0 background or 1 signal). Estimated using k-means.
#'  \item \code{mu} : The mean signal value in each state.
#'  \item \code{sigma} : The variances of the signal in each state. 
#' }
#' @export
#'
#' @examples
estimateInitialDomainsDistribution <- function(bw_file, seqlen, winsize=5000, stepsize=5000, blacklist=NULL, smooth=FALSE, w=3){
  require(rtracklayer)
  
  
  sizeDomain <- 0 ## The size in bp of regions in the 90th percentile enrichment
  sizeBground <- 0 ## The size in bp of regions in the background
  
  DomainScores <- c()
  BgrdScores <- c()
  
  for(ch in names(seqlen)){
    #cat("Processing " %+% red(ch) %+% "\n")
    region <- GRanges(ch, IRanges(1, seqlen[ch]))
    bins = slidingWindows(region,width = winsize,step = stepsize)[[1]]
    if(!is.null(blacklist)){
      bins = subsetByOverlaps(bins,blacklist,invert = T)
    }
    bins$score  = calcPeaksSignal(bins, bw_file)$meanScore
    
    q99 = quantile(bins$score, 0.99)
    m = max(bins$score)
    if(m/q99 > 3){
      bins$score[bins$score>q99] = q99
    }
    
    if(smooth){
      cat(paste0("Smoothing ",ch," scores"))
      bins$score = smoothSignal(bins$score,window = w)
    }
    
    set.seed(12345)
    km = kmeans(bins$score,centers = 2,nstart = 2,iter.max = 100)
    
    peaks_clus = which.max(km$centers)
    tmp1 = bins[which(km$cluster == peaks_clus)]
    tmp2 = bins[which(km$cluster != peaks_clus)]
    
    
    sizeDomain <- sizeDomain + sum(width(tmp1))
    sizeBground <- sizeBground + sum(width(tmp2))
    
    DomainScores <- c(DomainScores, tmp1$score)
    BgrdScores <- c(BgrdScores, tmp2$score)
    
  }

  
  
  P0 <- c(sizeDomain,sizeBground)/(sizeBground+sizeDomain)
  res <- list(P0=P0,
              mu= c(mean(DomainScores),mean(BgrdScores)),
              sigma=c(min(sd(DomainScores)^2,mean(DomainScores)/2) , min(sd(BgrdScores)^2,mean(BgrdScores)) ) )
  message("Initialization parameters")
  print(res)
  return(res)
}


#' Use the trained model to call domains (peaks)
#'
#' @param bw_file Path to the input bigWig file
#' @param hmmmodel An \LinkS4class(hmmspec) object specifiying the HMM model paramters.
#' @param genome The genome to use. It can be one of 'mm9','mm10' or 'hg19'.
#' @param chromosomes The chromosomes to call domains for.
#' @param prob_fout If you want to save the model probabilities as bigwig files set this value to a file path. Default: NULL.
#' @param winsize The bin size to use (in bp). Default: 5000.
#' @param stepsize The step size to use. Each time the window is slided by this value. Default: 5000.
#' @param smooth Do you want to smooth the signal before calling domains (peaks)? Default: FALSE.
#' @param w The number of bins used to smooth the signal. Default: 3. Ignored if smooth = FALSE. 
#' @param asDomains Do you want the function to only return the domains coordinates (TRUE) or you want the results at the bin level (FALSE)?
#'
#' @return
#' @export
#'
#' @examples
CallPeaksFromTrainedModel <- function(bw_file, hmmmodel,
                            genome = 'mm10',
                            chromosomes=NULL,
                            prob_fout =NULL,
                            winsize = 5e3,
                            stepsize = 5e3,
                            smooth = FALSE, w=3, 
                            asDomains = TRUE){
  
  refGenome = switch (genome,
                      mm9 = BSgenome.Mmusculus.UCSC.mm9::Mmusculus,
                      mm10 = BSgenome.Mmusculus.UCSC.mm10::Mmusculus,
                      hg19 = BSgenome.Mmusculus.UCSC.hg19::BSgenome.Mmusculus.UCSC.hg19
  )
  
  if(is.null(refGenome)){
    stop('the specified reference genome is not suported.')
  }
  
  if(is.null(chromosomes)){
    chromosomes = seqlevels(refGenome)
  }else{
    chromosomes = intersect(chromosomes, seqlevels(refGenome))
  }
  
  
  bw_scores <- list()
  
  # Call domains for each chromosome
  bins = GRangesList()
  pb <- txtProgressBar(min = 0, max = 4, style = 3)
  i=1
  for(ch in chromosomes){
    #cat(blue(paste("- chromosome",ch,"\n")))
    
    ## Binning chromosome
    region <- GRanges(ch, IRanges(1,seqlengths(refGenome)[ch]))
    bins[[ch]] = slidingWindows(region,width = winsize,step = stepsize)[[1]]
    
  }
  setTxtProgressBar(pb,1)
  bins = bins %>% GRangesList() %>% unlist()
  bins$scores  = calcPeaksSignal(bins, bw_file)$meanScore
  setTxtProgressBar(pb,2)
  ## trim scores if the max value is at least 5 folds larger than the 99th quantile
  m = max(bins$scores)
  q99 = quantile(bins$scores, 0.99)
  
  if(m/q99 > 5){
    bins$scores[bins$scores>q99] =q99
  }
  
  if(smooth){
    bins$scores = smoothSignal(bins$scores,window = w)
  }
  
  ## Predict enrichment
  train = list(x=bins$scores, N= length(bins$scores))
  y = predict(hmmmodel, train, method="smoothed")
  setTxtProgressBar(pb,3)
  # Define the states of the bins
  m1 = mean(bins$scores[y$s==1])
  m2 = mean(bins$scores[y$s==2])
  
  states = y$s
  if(m1>m2){
    states = ifelse(states==1 & y$p[,1]>=0,1,0)
    bins$score = 100 * y$p[,1]
  }else{
    states = ifelse(states==2 & y$p[,2]>=0,1,0)
    bins$score = 100 * y$p[,2]
  }
  
  bins$states= states
  if(!is.null(prob_fout)){
    cat(red("Saving domains probs\n"))
    seqlengths(bins) = seqlengths(refGenome)[chromosomes]
    export.bw(bins,con = prob_fout)
  }
  
  # Segment bins
  if(asDomains){
    tmp <- split(bins, as.character(states))
    tmp_merged = lapply(tmp,GenomicRanges::reduce)
    
    tmp_merged = lapply(c("0","1"), function(x) {tmp_merged[[x]]$score = as.numeric(x); tmp_merged[[x]]})
    
    tmp_merged = tmp_merged %>% GRangesList() %>% unlist()
    #tmp_merged$score = runValue(ttt)
    bw_scores = subset(tmp_merged,score==1)
    #bw_scores[[ch]] = tmp_merged
    
    bw_scores = unlist(GRangesList(bw_scores))
    setTxtProgressBar(pb,4)
    seqlengths(bw_scores) = seqlengths(refGenome)[chromosomes]
    return(bw_scores)
  }else{
    setTxtProgressBar(pb,4)
    seqlengths(bins) = seqlengths(refGenome)[chromosomes]
    return(bins)
  }
}


#' Calculates the mean bigWig signal per-bin
#'
#' @param windows A \LinkS4class{GRanges} object containing the The genomic coordinates of the bins.
#' @param bw Path to the bigWig file.
#'
#' @return
#' A \LinkS4class{GRanges} object with mean enrichement values in the column 'meanScore'.
#' @export
#'
#' @examples
calcPeaksSignal <- function(windows, bw){
  
  # check if windows have width > 1
  if( any(width(windows)==1) ){
    stop("provide 'windows' with widths greater than 1")
  }
  
  bwscores <- import.bw(bw, which= windows)
  covs = coverage(bwscores, weight=bwscores$score)
  covs = covs[seqlevels(windows)]
  windows <- GenomicRanges::binnedAverage(windows, covs, "meanScore")
  
  return(windows)
}


#' Smooths a numeric vector using a sliding window
#'
#' @param signal A numeric vector containing the values to smooth.
#' @param window The sliding window size we want to use. Should be an odd number. Default:3.
#'
#' @return
#' @export
#'
#' @examples
smoothSignal <- function(signal, window=3){
  
  smoothedVec = as.vector(stats::filter(signal, rep(1/(2*window+1), (2*window+1))))
  smoothedVec[1:window] = signal[1:window]
  smoothedVec[(length(signal)-window+1):length(signal)] = signal[(length(signal)-window+1):length(signal)]
  return(smoothedVec)
}



#' Main function used to estimate the HMM model parameters and call domains.
#'
#' @param bw_file Path to the input bigWig file
#' @param winsize The bin size to use (in bp). Default: 5000.
#' @param stepsize The step size to use. Each time the window is slided by this value. Default: 5000.
#' @param smooth Do you want to smooth the signal first then estimate paramters? Generally useful when the signal is sparse. Default: FALSE.
#' @param w The number of bins used to smooth the signal. Default: 3. Ignored if smooth = FALSE.
#' @param training.chrom The name of the chromosomes to use for training the model.
#' @param chromsToUse The name of the chromosomes to call the domains (peaks) for.
#' @param genome The genome to use. It can be 'mm9', 'mm10' or 'hg19'.
#' @param blacklist A \linkS4class{GRanges} object containing the coordinates of the blacklist regions to exclude. Default: NULL
#' @param saveProbs Whether to save the model posterior probabilities? Default: FALSE.
#' @param outDir Path to the output directory.
#'
#' @return
#' 
#' A \LinkS4class{GRanges} object containing the domains genomic locations.
#' @export
#'
#' @examples
CallDomains <- function(bw_file, 
                        winsize = 5000,
                        stepsize = 5000,
                        smooth = FALSE,
                        w=3,
                        training.chrom = glue("chr{1:4}"), 
                        chromsToUse = glue("chr{1:19}"),
                        genome="mm10", 
                        mm10_blacklist_gr=GRanges(),
                        saveProbs = FALSE,
                        plot.model=FALSE,
                        outDir = "Domains"                        
){
  
  dir.create(outDir,showWarnings = FALSE)
  refGenome = switch (genome,
                      mm9 = BSgenome.Mmusculus.UCSC.mm9::Mmusculus,
                      mm10 = BSgenome.Mmusculus.UCSC.mm10::Mmusculus,
                      hg19 = BSgenome.Mmusculus.UCSC.hg19::BSgenome.Mmusculus.UCSC.hg19
  )
  
  cat("\n")
  cat(bgGreen(glue("******** Processing {basename(bw_file)} ********")))
  cat(yellow("\n*) Training model\n"))
  hmmmodel = TrainHMMmodel(bw_file = bw_file,
                           seqlen = seqlengths(refGenome)[training.chrom],
                           blacklist = mm10_blacklist_gr,
                           plot.model = plot.model,
                           winsize = winsize,
                           step=stepsize,
                           smooth = smooth,
                           w=w
  )
  
  cat(yellow("*) Detecting domains:\n"))
  #cat(yellow(" ====================\n"))
  
  fprob = NULL
  if(saveProbs){
    fprob = basename(bw_file) %>% tools::file_path_sans_ext()
    fprob <- glue("{outDir}/{fprob}_{round(winsize/1e3,2)}_probs.bw")
  }

  Domains.gr = suppressWarnings(CallPeaksFromTrainedModel(bw_file, 
                                                hmmmodel, 
                                                genome = 'mm10',
                                                chromosomes=chromsToUse,
                                                prob_fout = fprob,
                                                smooth = smooth,
                                                w = w))
  
  cat(yellow("\nGenerating bed file\n"))
  cat(yellow(" ====================\n"))
  fout = basename(bw_file) %>% tools::file_path_sans_ext()
  fout <- glue("{outDir}/{fout}_{round(winsize/1e3,2)}_domains.bed")
  
  export.bed(Domains.gr,con = fout)
  cat(green(paste0("Saved in :", fout,"\n")))
  cat("\n")
  return(Domains.gr)
}
