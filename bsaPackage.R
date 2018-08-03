#This file is used to smooth the data and plot them. Code is modified version of the QTLSeqR package




# Function to set comparisons between two samples and output the table in a way that the QTLseq package can use.

#' Uses two element of a list to create a table that can be used by the package QTLseq.
#' 
#' @param highBulk String. The name of the sample (and of the element in the SNPset list) corresponding to the High Bulk.
#' 
#' @return a data.frame with the format that can be understood by QTLseq. 

picker=function(sequenced){
  # if(!lowBulk %in% names(SNPset) | !highBulk %in% names(SNPset)){
  #   stop("Confirm if the bulks' names are present in the SNP set (vcf file) provided.")
  # }
  metaTable=data.frame(CHROM=sequenced$CHR, POS=sequenced$POS, REF=NA, 
                       ALT=NA, 
                       AD_REF.LOW=sequenced$ref1,
                       AD_ALT.LOW=sequenced$alt1,
                       DP.LOW=sequenced$depth1,
                       GQ.LOW=NA ,
                       PL.LOW=NA,
                       SNPindex.LOW=sequenced$SNPindL,
                       AD_REF.HIGH=sequenced$ref2,
                       AD_ALT.HIGH=sequenced$alt2,
                       DP.HIGH=sequenced$depth2,
                       GQ.HIGH=NA ,
                       PL.HIGH=NA,
                       SNPindex.HIGH= sequenced$SNPindH, #AD_ALT.HIGH/DP.HIGH,
                       REF_FRQ=((sequenced$ref1 + sequenced$ref2)/(sequenced$depth1 + sequenced$depth2))
  )
  metaTable=dplyr::mutate(metaTable, deltaSNP=SNPindex.HIGH-SNPindex.LOW)
  return(metaTable)
}

SNPset<- picker(sequenced)


#' @param SNPset The data frame imported by \code{ImportFromGATK} 
#' @param windowSize the window size (in base pairs) bracketing each SNP for which to calculate the statitics.
#' @param popStruc the population structure. Defaults to "F2" and assumes "RIL" otherwise
#' @param bulkSize non-negative integer vector. The number of individuals in
#'   each simulated bulk. Can be of length 1, then both bulks are set to the
#'   same size. Assumes the first value in the vector is the simulated high
#'   bulk.
#' @param depth integer. A read depth for which to replicate SNP-index calls.
#' @param replications integer. The number of bootstrap replications.
#' @param filter numeric. An optional minimum SNP-index filter
#' @param intervals numeric vector. Confidence intervals supplied as two-sided percentiles. i.e. If intervals = '95' will return the two sided 95\% confidence interval, 2.5\% on each side. 
#'
#' @return A SNPset data frame with delta SNP-index thresholds corrisponding to the 
#' requested confidence intervals matching the tricube smoothed depth at each SNP.
#' @export runQTLseqAnalysis
#'
#' @examples df_filt <- runQTLseqAnalysis(
#' SNPset = df_filt,
#' bulkSize = c(25, 35)
#' windowSize = 1e6,
#' popStruc = "F2",
#' replications = 10000,
#' intervals = c(95, 99)
#' )
#' 

runQTLseqAnalysis2 <- function(SNPset,
                              windowSize = 25000,
                              
                              bulkSize,
                              depth = NULL
                              
                             
                              ) {
  
  message("Counting SNPs in each window...")
  SNPset <- SNPset %>%
    dplyr::group_by(CHROM) %>%
    dplyr::mutate(nSNPs = countSNPs_cpp(POS = POS, windowSize = windowSize))
  
  message("Calculating tricube smoothed delta SNP index...")
  SNPset <- SNPset %>%
    dplyr::mutate(tricubeDeltaSNP = tricubeStat(POS = POS, Stat = deltaSNP, windowSize))
  
  # #convert intervals to quantiles
  # if (all(intervals >= 1)) {
  #   message(
  #     "Returning the following two sided confidence intervals: ",
  #     paste(intervals, collapse = ", ")
  #   )
  #   quantiles <- (100 - intervals) / 200
  # } else {
  #   stop(
  #     "Convidence intervals ('intervals' paramater) should be supplied as two-sided percentiles. i.e. If intervals = '95' will return the two sided 95% confidence interval, 2.5% on each side."
  #   )
  # }
  # 
  #calculate min depth per snp between bulks
  SNPset <-
    SNPset %>%
    dplyr::mutate(minDP = pmin(DP.LOW, DP.HIGH))
  
  SNPset <-
    SNPset %>%
    dplyr::group_by(CHROM) %>%
    dplyr::mutate(tricubeDP = floor(tricubeStat(POS, minDP, windowSize = windowSize)))
  
  if (is.null(depth)) {
    message(
      "Variable 'depth' not defined, using min and max depth from data: ",
      min(SNPset$minDP),
      "-",
      max(SNPset$minDP)
    )
    depth <- min(SNPset$minDP):max(SNPset$minDP)
  }
  # 
  # #simualte confidence intervals
  # CI <-
  #   simulateConfInt(
  #     popStruc = popStruc,
  #     bulkSize = bulkSize,
  #     depth = depth,
  #     replications = replications,
  #     filter = filter,
  #     intervals = quantiles
  #   )
  # 
  # 
  # #match name of column for easier joining of repeat columns
  # names(CI)[1] <- "tricubeDP"
  # 
  #use join as a quick way to match min depth to matching conf intervals.
  # SNPset <-
  #   dplyr::left_join(x = SNPset,
  #                    y = CI #, commented out becuase of above change. need to remove eventually
  #                    # by = c("tricubeDP" = "depth"))
  #   )
  return(as.data.frame(SNPset))
  
}

tricubeStat <- function(POS, Stat, windowSize = 2e6)
{
  if (windowSize <= 0)
    stop("A positive smoothing window is required")
  stats::predict(locfit::locfit(Stat ~ locfit::lp(POS, h = windowSize, deg = 0)), POS)
}

plotQTLStats2 <-
  function(SNPset,
           subset = NULL,
           var = "nSNPs",
           scaleChroms = TRUE,
           line = TRUE,
           plotThreshold = FALSE,
           plotIntervals = FALSE,
           q = 0.05,
           ...) {
    
    #get fdr threshold by ordering snps by pval then getting the last pval
    #with a qval < q
    
    if (!all(subset %in% unique(SNPset$CHROM))) {
      whichnot <-
        paste(subset[base::which(!subset %in% unique(SNPset$CHROM))], collapse = ', ')
      stop(paste0("The following are not true chromosome names: ", whichnot))
    }
    
    if (!var %in% c("nSNPs", "deltaSNP", "Gprime", "negLog10Pval"))
      stop(
        "Please choose one of the following variables to plot: \"nSNPs\", \"deltaSNP\", \"Gprime\", \"negLog10Pval\""
      )
    
    #don't plot threshold lines in deltaSNPprime or number of SNPs as they are not relevant
    if ((plotThreshold == TRUE &
         var == "deltaSNP") |
        (plotThreshold == TRUE & var == "nSNPs")) {
      message("FDR threshold is not plotted in deltaSNP or nSNPs plots")
      plotThreshold <- FALSE
    }
    #if you need to plot threshold get the FDR, but check if there are any values that pass fdr
    
    GprimeT <- 0
    logFdrT <- 0
    
    if (plotThreshold == TRUE) {
      fdrT <- getFDRThreshold(SNPset$pvalue, alpha = q)
      
      if (is.na(fdrT)) {
        warning("The q threshold is too low. No threshold line will be drawn")
        plotThreshold <- FALSE
        
      } else {
        logFdrT <- -log10(fdrT)
        GprimeT <- SNPset[which(SNPset$pvalue == fdrT), "Gprime"]
      }
    }
    
    SNPset <-
      if (is.null(subset)) {
        SNPset
      } else {
        SNPset[SNPset$CHROM %in% subset,]
      }
    
    p <- ggplot2::ggplot(data = SNPset) +
      ggplot2::scale_x_continuous(breaks = seq(from = 0,to = max(SNPset$POS), by = 10^(floor(log10(max(SNPset$POS))))), labels = format_genomic(), name = "Genomic Position (Mb)") +
      ggplot2::theme(plot.margin = ggplot2::margin(
        b = 10,
        l = 20,
        r = 20,
        unit = "pt"
      ))
    
    if (var == "Gprime") {
      threshold <- GprimeT
      p <- p + ggplot2::ylab("G' value")
    }
    
    if (var == "negLog10Pval") {
      threshold <- logFdrT
      p <-
        p + ggplot2::ylab(expression("-" * log[10] * '(p-value)'))
    }
    
    if (var == "nSNPs") {
      p <- p + ggplot2::ylab("Number of SNPs in window")
    }
    
    if (var == "deltaSNP") {
      var <- "tricubeDeltaSNP"
      p <-
        p + ggplot2::ylab(expression(Delta * '(SNP-index)')) +
        ggplot2::ylim(-0.55, 0.55) +
        ggplot2::geom_hline(yintercept = 0,
                            color = "black",
                            alpha = 0.4)
      if (plotIntervals == TRUE) {
        
        ints_df <-
          dplyr::select(SNPset, CHROM, POS, dplyr::matches("CI_")) %>% tidyr::gather(key = "Interval", value = "value",-CHROM,-POS)
        
        p <- p + ggplot2::geom_line(data = ints_df, ggplot2::aes(x = POS, y = value, color = Interval)) +
          ggplot2::geom_line(data = ints_df, ggplot2::aes(
            x = POS,
            y = -value,
            color = Interval
          ))
      }
    }
    
    if (line) {
      p <-
        p + ggplot2::geom_line(ggplot2::aes_string(x = "POS", y = var), ...)
    }
    
    if (!line) {
      p <-
        p + ggplot2::geom_point(ggplot2::aes_string(x = "POS", y = var), ...)
    }
    
    if (plotThreshold == TRUE)
      p <-
      p + ggplot2::geom_hline(
        ggplot2::aes_string(yintercept = "threshold"),
        color = "red",
        size = 1,
        alpha = 0.4
      )
    
    if (scaleChroms == TRUE) {
      p <- p + ggplot2::facet_wrap(~ CHROM, scales = "free_x")
    } else {
      p <- p + ggplot2::facet_wrap(~ CHROM, scales = "free_x")    
    }
    
    p
  }

format_genomic <- function(...) {
  # Format a vector of numeric values according
  # to the International System of Units.
  # http://en.wikipedia.org/wiki/SI_prefix
  #
  # Based on code by Ben Tupper
  # https://stat.ethz.ch/pipermail/r-help/2012-January/299804.html
  # Args:
  #   ...: Args passed to format()
  #
  # Returns:
  #   A function to format a vector of strings using
  #   SI prefix notation
  #
  
  function(x) {
    limits <- c(1e0,   1e3, 1e6)
    #prefix <- c("","Kb","Mb")
    
    # Vector with array indices according to position in intervals
    i <- findInterval(abs(x), limits)
    
    # Set prefix to " " for very small values < 1e-24
    i <- ifelse(i==0, which(limits == 1e0), i)
    
    paste(format(round(x/limits[i], 1),
                 trim=TRUE, scientific=FALSE, ...)
          #  ,prefix[i]
    )
  }
}


SNPset <- runQTLseqAnalysis2(SNPset, bulkSize = 10 , windowSize = 25000 )


plotQTLStats2(SNPset, var = "deltaSNP", scaleChroms = TRUE)

