# This file contains the functions for the parallel peak filling step in the XCMS workflow
#' @importFrom methods setMethod validObject
.xmse_process_history <- utils::getFromNamespace(".xmse_process_history", "xcms")
.features_ms_region <- utils::getFromNamespace(".features_ms_region", "xcms")
.history2fill_fun <- utils::getFromNamespace(".history2fill_fun", "xcms")
.xmse_integrate_chrom_peaks <- utils::getFromNamespace(".xmse_integrate_chrom_peaks", "xcms")
.subset_xcms_experiment <- utils::getFromNamespace(".subset_xcms_experiment", "xcms")
.PROCSTEP.PEAK.DETECTION <- utils::getFromNamespace(".PROCSTEP.PEAK.DETECTION", "xcms")
.PROCSTEP.PEAK.FILLING <- utils::getFromNamespace(".PROCSTEP.PEAK.FILLING", "xcms")
.chromPeaks <- utils::getFromNamespace(".chromPeaks", "xcms")
.featureIDs <- utils::getFromNamespace(".featureIDs", "xcms")

#' Preparation for the computational intensive peak filling step
#' @description cPrepFillChromPeaks, cGetChromPeakData, and cFillChromPeaks are small constructs of fillChromPeaks.
#' @description By running the functions as suggested in the vignette, we minimize the RAM usage, thus making the workflow scalable to tens of thousands of files.
#' @param input path to XcmsExperiment as returned by the peak alignment step
#' @param output path(s) of the output file(s). This file goes into cFillChromPeaksStep2()
#' @param param the object returned by xcms::FillChromPeaksParam()
#' @param interval numeric vector of length 2 with start and end sample indices
#' @param settings settings list from settings.yaml
#' @import xcms
#' @importFrom MSnbase fileNames requiredFvarLabels
#' @importFrom readr read_rds write_rds
#' @importFrom purrr imap
#' @return a list object with all the necessary objects for the cGetChromPeakData and cFillChromPeaks functions
#' @export
cFillChromPeaksStep1 <- function(input, output, param, interval, settings){

  msLevel = 1L
  object <- readr::read_rds(input)

  feature_ids <- rownames(featureDefinitions(object, msLevel = msLevel))
  fr <- .features_ms_region(object, mzmin = param@mzmin,
                            mzmax = param@mzmax, rtmin = param@rtmin,
                            rtmax = param@rtmax, features = feature_ids)

  fr <- cbind(fr, mzmed = featureDefinitions(object, msLevel = msLevel)$mzmed)
  fvals <- featureValues(object, value = "index", msLevel = msLevel)
  pal <- lapply(seq_len(ncol(fvals)), function(i) {fr[is.na(fvals[, i]), , drop = FALSE]})
  names(pal) <- seq_along(pal)

  ph <- .xmse_process_history(object, .PROCSTEP.PEAK.DETECTION, msLevel = msLevel)
  fill_fun <- .history2fill_fun(ph)
  mzf <- "wMean"
  if (length(ph) && inherits(ph[[1L]], "XProcessHistory")) {
    prm <- ph[[1L]]@param
    if (any(slotNames(prm) == "mzCenterFun"))
      mzf <- prm@mzCenterFun
  } else
    prm <- MatchedFilterParam()
  mzf <- paste0("mzCenter.", gsub("mzCenter.", "", mzf, fixed = TRUE))

  seq_along(object)[interval[1]:interval[2]] |> purrr::imap(~{
    if (is.na(.x)) return() #  | file.exists(output[.y])

    run_package <- list(
      "object" = .subset_xcms_experiment(
        object,
        i = .x,
        keepAdjustedRtime = TRUE,
        ignoreHistory = TRUE),
      "pal" = pal[.x],
      "intFun" = fill_fun,
      "mzCenterFun" = mzf,
      "param" = prm
    )
    readr::write_rds(run_package, file = output[.y])
  })
  return()
}



singleSampleFillChromPeaks <- function(object, pal, intFun, mzCenterFun, param, output){
  res <-
    .xmse_integrate_chrom_peaks(
      x = object,
      pal = pal,
      intFun = intFun,
      mzCenterFun = mzCenterFun,
      param = param)
  
  readr::write_rds(res, output)
  return()
}

#' Subset XcmsExperiment
#'
#' @name [
#' @aliases [,XcmsExperiment,ANY,ANY,ANY-method
#' @param x XcmsExperiment object
#' @param i row indices
#' @param j column indices
#' @param ... additional arguments
#' @param drop logical, whether to drop dimensions
#' @return subsetted XcmsExperiment
#' @keywords internal
setMethod("[", "XcmsExperiment", function(x, i, j, ..., drop = TRUE) {
  callNextMethod()
})

#' Intensive peak extraction step performed file-wise (parallel in cluster workflows)
#'
#' @param input path to the prepared data file from cFillChromPeaksStep1
#' @param output the output file path for storage of extracted peaks
#' @import xcms
#' @importFrom readr read_rds
#' @importFrom readr write_rds
#' @return NULL (writes output to file)
#' @export
cFillChromPeaksStep2 <- function(input, output){
  prepared_data <- readr::read_rds(input)
  prepared_data$output <- output
  do.call(singleSampleFillChromPeaks, prepared_data)
  return()
  }

.feature_values <- function(pks, fts, method, value = "into",
                            intensity = "into", colnames,
                            missing = NA) {
  ftIdx <- fts$peakidx
  ## Match columns
  idx_rt <- match("rt", colnames(pks))
  idx_int <- match(intensity, colnames(pks))
  idx_samp <- match("sample", colnames(pks))
  vals <- matrix(nrow = length(ftIdx), ncol = length(colnames))
  nSamples <- seq_along(colnames)
  if (method == "sum") {
    for (i in seq_along(ftIdx)) {
      cur_pks <- pks[ftIdx[[i]], c(value, "sample"), drop=FALSE]
      int_sum <- split(as.numeric(cur_pks[, value]),
                       as.factor(as.integer(cur_pks[, "sample"])))
      vals[i, as.numeric(names(int_sum))] <-
        unlist(lapply(int_sum, base::sum), use.names = FALSE)
    }
  } else {
    if (method == "medret") {
      medret <- fts$rtmed
      for (i in seq_along(ftIdx)) {
        gidx <- ftIdx[[i]][
          base::order(base::abs(as.numeric(pks[ftIdx[[i]],
                                    idx_rt]) - medret[i]))]
        vals[i, ] <- gidx[
          base::match(nSamples, as.numeric(pks[gidx, idx_samp]))]
      }
    }
    if (method == "maxint") {
      for (i in seq_along(ftIdx)) {
        gidx <- ftIdx[[i]][
          base::order(as.numeric(pks[ftIdx[[i]], idx_int]),
                      decreasing = TRUE)]
        vals[i, ] <- gidx[base::match(nSamples,
                                      as.numeric(pks[gidx, idx_samp]))]
      }
    }
    if (value != "index") {
      if (!any(colnames(pks) == value))
        stop("Column '", value, "' not present in the ",
             "chromatographic peaks matrix!")
      vals <- as.numeric(pks[vals, value])
      dim(vals) <- c(length(ftIdx), length(nSamples))
    }
  }
  if (value != "index") {
    if (is.numeric(missing)) {
      vals[is.na(vals)] <- missing
    }
    if (!is.na(missing) & missing == "rowmin_half") {
      for (i in seq_len(nrow(vals))) {
        nas <- is.na(vals[i, ])
        if (any(nas))
          vals[i, nas] <- min(vals[i, ], na.rm = TRUE) / 2
      }
    }
  }
  colnames(vals) <- colnames
  rownames(vals) <- rownames(fts)
  vals
}

#' fill chrom peak function 3
#'
#' @param inputFromStep1 Output of cFillChromPeaksStep1() function
#' @param inputFromStep2 Output of cFillChromPeaksStep2() function
#' @param output the output which contains the result of the preprocessing
#' @param param FillChromPeaksParam object
#' @param settings settings list from settings.yaml
#'
#' @return a preprocessed xcms object
#' @import xcms
#' @importFrom MSnbase fileNames
#' @importFrom readr read_rds read_csv
#' @importFrom S4Vectors extractROWS
#' @importFrom purrr map walk2
#' @importFrom progressr with_progress progressor
#' @export
cFillChromPeaksStep3 <- function(inputFromStep1, inputFromStep2, output, param,
                                 settings) {

  object <- readr::read_rds(inputFromStep1)
  msLevel <- 1
  run_scheduler <- readr::read_csv(paste0(settings$general$output_path, "run_schedule_", settings$general$run_id, ".csv"))

  res <- run_scheduler$integrated2 |> purrr::map(readr::read_rds, .progress = TRUE)
  res <- do.call(rbind, res)

  i_res <- seq((nrow(.chromPeaks(object)) + 1L), length.out = nrow(res))
  i_res <- split(i_res, rownames(res))
  i_ft <- match(names(i_res), rownames(featureDefinitions(object)))

  with_progress({
    p <- progressor(along = i_res)

    walk2(
      .x = i_res,
      .y = i_ft,
      .f = function(idx, ft) {
        if (is.na(ft)) {
          warning("Feature ID not found in featureDefinitions")
          p()
          return(NULL)
        }
        object@featureDefinitions$peakidx[[ft]] <<-
          sort(c(object@featureDefinitions$peakidx[[ft]], idx))
        p()
      }
    )
  })

  nr <- nrow(res)
  maxi <- max(as.integer(sub("CP", "", rownames(.chromPeaks(object)))))
  rownames(res) <- .featureIDs(nr, "CP", maxi + 1)
  cpd <- data.frame(ms_level = rep(as.integer(msLevel), nr),
                    is_filled = rep(TRUE, nr))
  rownames(cpd) <- rownames(res)
  object@chromPeaks <- rbind(object@chromPeaks, res)
  object@chromPeakData <- MsCoreUtils::rbindFill(object@chromPeakData, cpd)

  ph <- xcms:::XProcessHistory(param = param,
                               date. = date(),
                               type. = .PROCSTEP.PEAK.FILLING,
                               fileIndex = seq_along(object),
                               msLevel = msLevel)
  object <- xcms:::addProcessHistory(object, ph)
  validObject(object)

  return(object)
}
