
.getPeakGroupsRtMatrix <- utils::getFromNamespace(".getPeakGroupsRtMatrix", "xcms")
.peakIndex <- utils::getFromNamespace(".peakIndex", "xcms")
.getPeakGroupsRtMatrix <- utils::getFromNamespace(".getPeakGroupsRtMatrix", "xcms")
na.flatfill <- utils::getFromNamespace("na.flatfill", "xcms")
adjustRtimeSubset <- utils::getFromNamespace("adjustRtimeSubset", "xcms")
processHistory <- utils::getFromNamespace("processHistory", "xcms")
.has_chrom_peak_data <- utils::getFromNamespace(".has_chrom_peak_data", "xcms")
.PROCSTEP.PEAK.GROUPING <- utils::getFromNamespace(".PROCSTEP.PEAK.GROUPING", "xcms")
.PROCSTEP.PEAK.DETECTION <- utils::getFromNamespace(".PROCSTEP.PEAK.DETECTION", "xcms")
.PROCSTEP.RTIME.CORRECTION <- utils::getFromNamespace(".PROCSTEP.RTIME.CORRECTION", "xcms")



#' Auxillary function for adjust R time function
#'
#' @param object ms object
#' @param param xcms::peakgroupsparam object
#'
#' @return pkGrp object
adjustRtimePeakGroups <- function(object, param = PeakGroupsParam()) {

  msLevel = 1L
  if (!is(object, "XCMSnExp"))
    stop("'object' has to be an 'XCMSnExp' object.")
  if (!hasFeatures(object))
    stop("No features present. Please run 'groupChromPeaks' first.")
  if (hasAdjustedRtime(object))
    warning("Alignment/retention time correction was already performed, ",
            "returning a matrix with adjusted retention times.")

  subs <- integer(0)
  if (!length(subs)) {subs <- seq_along(fileNames(object))}

  nSamples <- length(subs)
  pkGrp <- .getPeakGroupsRtMatrix(
    peaks = chromPeaks(object, msLevel = msLevel),
    peakIndex = .peakIndex(featureDefinitions(object)), # This line is changed. We do not update feature definitions.
    sampleIndex = subs,
    missingSample = nSamples - (nSamples * minFraction(param)),
    extraPeaks = extraPeaks(param)
  )
  colnames(pkGrp) <- basename(fileNames(object))[subs]
  return(pkGrp)
}



#' Title
#'
#' @param peaks The identified peaks and their positions
#' @param peakIndex index of peak
#' @param rtime retention times of samples
#' @param minFraction the minimum number of samples a feature has to exist in to qualify as an "aligning feature"
#' @param extraPeaks How many extra peaks are allowed in the area of alignment (can be high for large datasets with high drift)
#' @param smooth smoothing of alignment
#' @param span how many features should the non-linear alignment model look at simultaniously? close to 0 models highly non-linear drift and close to 1 models linear drift
#' @param family kind of model
#' @param peakGroupsMatrix the groups (used for minFraction)
#' @param subset The code is rewritten such that the subset option is no longer available. All samples will be aligned.
#' @param subsetAdjust see "subset"
#'
#' @return adjusted rtimes
do_adjustRtime_peakGroups <-
  function(peaks, peakIndex, rtime, minFraction = 0.9, extraPeaks = 1,
           smooth = c("loess", "linear"), span = 0.2,
           family = c("gaussian", "symmetric"),
           peakGroupsMatrix = matrix(ncol = 0, nrow = 0),
           subset = integer(), subsetAdjust = c("average", "previous")){

    ## Translate minFraction to number of allowed missing samples.
    nSamples      <- length(subset)
    missingSample <- nSamples - (nSamples * minFraction)
    if (nrow(peakGroupsMatrix)) {
      rt <- peakGroupsMatrix
    } else
      rt <- .getPeakGroupsRtMatrix(peaks, peakIndex, subset,
                                   missingSample)

    message("Performing retention time correction using ", nrow(rt),
            " peak groups.")

    ## Calculate the deviation of each peak group in each sample from its
    ## median
    rtdev <- rt - apply(rt, 1, stats::median, na.rm = TRUE)

    if (smooth == "loess") {
      mingroups <- min(colSums(!is.na(rt)))
      if (mingroups < 4) {
        smooth <- "linear"
        warning("Too few peak groups for 'loess', reverting to linear",
                " method")
      } else if (mingroups * span < 4) {
        span <- 4 / mingroups
        warning("Span too small for 'loess' and the available number of ",
                "peak groups, resetting to ", round(span, 2))
      }
    }

    rtdevsmo <- vector("list", nSamples)

    ## Code for checking to see if retention time correction is overcorrecting
    rtdevrange <- range(rtdev, na.rm = TRUE)
    warn.overcorrect <- FALSE
    warn.tweak.rt <- FALSE

    rtime_adj <- rtime
    ## Adjust samples in subset.
    for (i in seq_along(subset)) {
      i_all <- subset[i]              # Index of sample in whole dataset.
      pts <- stats::na.omit(data.frame(rt = rt[, i], rtdev = rtdev[, i]))

      ## order the data.frame such that rt and rtdev are increasingly ordered.
      pk_idx <- order(pts$rt, pts$rtdev)
      pts <- pts[pk_idx, ]
      if (smooth == "loess") {
        lo <- suppressWarnings(stats::loess(rtdev ~ rt, pts, span = span,
                                     degree = 1, family = family))

        rtdevsmo[[i]] <- na.flatfill(
          stats::predict(lo, data.frame(rt = rtime[[i_all]])))
        ## Remove singularities from the loess function
        rtdevsmo[[i]][abs(rtdevsmo[[i]]) >
                        stats::quantile(abs(rtdevsmo[[i]]), 0.9,
                                 na.rm = TRUE) * 2] <- NA
        if (length(naidx <- which(is.na(rtdevsmo[[i]])))){
          rtdevsmo[[i]][naidx] <- suppressWarnings(
            stats::approx(stats::na.omit(data.frame(rtime[[i_all]], rtdevsmo[[i]])),
                   xout = rtime[[i_all]][naidx], rule = 2)$y
          )
          }

        ## Check if there are adjusted retention times that are not ordered
        ## increasingly. If there are, search for each first unordered rt
        ## the next rt that is larger and linearly interpolate the values
        ## in between (see issue #146 for an illustration).
        while (length(decidx <- which(diff(rtime[[i_all]] - rtdevsmo[[i]]) < 0))) {
          warn.tweak.rt <- TRUE  ## Warn that we had to tweak the rts.
          rtadj <- rtime[[i_all]] - rtdevsmo[[i]]
          rtadj_start <- rtadj[decidx[1]] ## start interpolating from here
          ## Define the
          next_larger <- which(rtadj > rtadj[decidx[1]])
          if (length(next_larger) == 0) {
            ## Fix if there is no larger adjusted rt up to the end.
            next_larger <- length(rtadj) + 1
            rtadj_end <- rtadj_start
          } else {
            next_larger <- min(next_larger)
            rtadj_end <- rtadj[next_larger]
          }
          ## linearly interpolate the values in between.
          adj_idxs <- (decidx[1] + 1):(next_larger - 1)
          incr <- (rtadj_end - rtadj_start) / length(adj_idxs)
          rtdevsmo[[i]][adj_idxs] <- rtime[[i_all]][adj_idxs] -
            (rtadj_start + (1:length(adj_idxs)) * incr)
        }

        rtdevsmorange <- range(rtdevsmo[[i]])
        if (any(rtdevsmorange / rtdevrange > 2))
          warn.overcorrect <- TRUE
      } else {
        if (nrow(pts) < 2) {
          stop("Not enough peak groups even for linear smoothing ",
               "available!")
        }
        ## Use lm instead?
        fit <- stats::lsfit(pts$rt, pts$rtdev)
        rtdevsmo[[i]] <- rtime[[i_all]] * fit$coef[2] + fit$coef[1]
        ptsrange <- range(pts$rt)
        minidx <- rtime[[i_all]] < ptsrange[1]
        maxidx <- rtime[[i_all]] > ptsrange[2]
        rtdevsmo[[i]][minidx] <- rtdevsmo[[i]][head(which(!minidx), n = 1)]
        rtdevsmo[[i]][maxidx] <- rtdevsmo[[i]][tail(which(!maxidx), n = 1)]
      }
      ## Finally applying the correction
      rtime_adj[[i_all]] <- rtime[[i_all]] - rtdevsmo[[i]]
    }
    ## Adjust the remaining samples.
    rtime_adj <- adjustRtimeSubset(rtime, rtime_adj, subset = subset, method = subsetAdjust)
    return(rtime_adj)
  }

#' Adjust retention time
#'
#' @param object ms object
#' @param param xcms::PeakGroupsParam
#'
#' @return aligned ms object
#' @export
adjustRtime <-
  function(object, param) {
            if (hasChromPeaks(object) & !.has_chrom_peak_data(object))
              object <- updateObject(object)
            msLevel <- 1L
            startDate <- date()
            ## If param does contain a peakGroupsMatrix extract that one,
            ## otherwise generate it.
            if (nrow(peakGroupsMatrix(param))){
              pkGrpMat <- peakGroupsMatrix(param)
              } else {
              pkGrpMat <- adjustRtimePeakGroups(object, param = param)
              }
            message("Data must only contain one(!) MS level")
            res <- do_adjustRtime_peakGroups(
              peaks = chromPeaks(object, msLevel = msLevel),
              # OBS: removed .update_feature_definitions, since we only have one MS level!
              peakIndex = featureDefinitions(object)$peakidx,
              rtime = rtime(object, bySample = TRUE),
              minFraction = minFraction(param),
              extraPeaks = extraPeaks(param),
              smooth = "loess",
              span = 0.8,
              family = "gaussian",
              peakGroupsMatrix = pkGrpMat,
              # OBS: subset is always the full dataset
              subset = 1:length(fileNames(object)),
              subsetAdjust = subsetAdjust(param)
            )
            ## Add the pkGrpMat that's being used to the param object.
            peakGroupsMatrix(param) <- pkGrpMat
            ## Dropping the peak groups but don't remove its process history
            ## step.
            ph <- processHistory(object, type = .PROCSTEP.PEAK.GROUPING)
            object <- dropFeatureDefinitions(object)
            ## Add the results. adjustedRtime<- should also fix the retention
            ## times for the peaks! Want to keep also the latest alignment
            ## information
            adjustedRtime(object) <- res
            if (length(ph)) {
              object <- addProcessHistory(object, ph[[length(ph)]])
            }
            ## Add the process history step, get the msLevel from the peak
            ## detection step.
            ph <- processHistory(object, type = .PROCSTEP.PEAK.DETECTION)
            xph <- XProcessHistory(param = param, date. = startDate,
                                   type. = .PROCSTEP.RTIME.CORRECTION,
                                   fileIndex = 1:length(fileNames(object)),
                                   msLevel = msLevel)
            object <- addProcessHistory(object, xph)
            validObject(object)
            return(object)
  }
