

.hasFilledPeaks   <- utils::getFromNamespace(".hasFilledPeaks", "xcms")
.param2string     <- utils::getFromNamespace(".param2string", "xcms")
.getChromPeakData <- utils::getFromNamespace(".getChromPeakData", "xcms")
.copy_env         <- utils::getFromNamespace(".copy_env", "xcms")
XProcessHistory   <- utils::getFromNamespace("XProcessHistory", "xcms")
addProcessHistory <- utils::getFromNamespace("addProcessHistory", "xcms")
.get_closest_index <- utils::getFromNamespace(".get_closest_index", "xcms")

#' Preparation for the computational intensive peak filling step
#' @description cPrepFillChromPeaks, cGetChromPeakData, and cFillChromPeaks are small constructs of fillChromPeaks.
#' @description By running the functions as suggested in the vignette, we minimize the RAM usage, thus making the workflow scalable to tens of thousands of files.
#' @param input path to XCMSset as returned by the peak alignment step
#' @param fcp the object returned by xcms::FillChromPeaksParam()
#' @param output path of the output file. This file goes into cFillChromPeaksStep2()
#' @import xcms
#' @importFrom stats median
#' @importFrom MSnbase fileNames
#' @importFrom MSnbase requiredFvarLabels
#' @importFrom readr read_rds
#' @importFrom readr write_rds
#' @return a list object with all the necessary objects for the cGetChromPeakData and cFillChromPeaks functions
#' @export
cFillChromPeaksStep1 <- function(input, output, fcp){
  msLevel = 1L
  object <- readr::read_rds(input)

  if (!xcms::hasFeatures(object)){stop("No feature definitions for MS level ", msLevel," present. Please run 'groupChromPeaks' first.")}
  if (.hasFilledPeaks(object)) {message("Filled peaks already present, adding still missing", " peaks.")}

  expandMz  <- xcms::expandMz(fcp)
  expandRt  <- xcms::expandRt(fcp)
  fixedMz   <- xcms::fixedMz(fcp)
  fixedRt   <- xcms::fixedRt(fcp)
  ppm       <- xcms::ppm(fcp)

  ## Define or extend the peak area from which the signal should be
  ## extracted.
  fdef       <- xcms::featureDefinitions(object)
  tmp_pks    <- xcms::chromPeaks(object)[, c("rtmin", "rtmax", "mzmin","mzmax")]
  aggFunLow  <- stats::median
  aggFunHigh <- stats::median

  pkArea <- do.call(
    rbind,
    lapply(
      fdef$peakidx,
      function(z) {
        pa <- c(
          aggFunLow(tmp_pks[z, 1]), aggFunHigh(tmp_pks[z, 2]), # Retention times
          aggFunLow(tmp_pks[z, 3]), aggFunHigh(tmp_pks[z, 4])  # m/z ratios
        )
        ## Check if we have to apply ppm replacement:
        if (ppm != 0) {
          mzmean <- mean(pa[3:4])
          tittle <- mzmean * (ppm / 2) / 1E6
          if ((pa[4] - pa[3]) < (tittle * 2)) {
            pa[3] <- mzmean - tittle
            pa[4] <- mzmean + tittle
          }
        }
        ## Expand it.
        if (expandRt != 0) {
          diffRt <- (pa[2] - pa[1]) * expandRt / 2
          pa[1] <- pa[1] - diffRt
          pa[2] <- pa[2] + diffRt
        }
        if (expandMz != 0) {
          diffMz <- (pa[4] - pa[3]) * expandMz / 2
          pa[3] <- pa[3] - diffMz
          pa[4] <- pa[4] + diffMz
        }
        if (fixedMz != 0) {
          pa[3] <- pa[3] - fixedMz
          pa[4] <- pa[4] + fixedMz
        }
        if (fixedRt != 0) {
          pa[1] <- pa[1] - fixedRt
          pa[2] <- pa[2] + fixedRt
        }
        pa
      }
    ))
  rm(tmp_pks)

  ## Add mzmed column - needed for MSW peak filling.
  colnames(pkArea) <- c("rtmin", "rtmax", "mzmin", "mzmax")
  pkArea           <- cbind(group_idx = 1:nrow(pkArea),
                            pkArea,
                            mzmed = as.numeric(fdef$mzmed))
  pkGrpVal         <- xcms::featureValues(object, value = "index")

  ## Check if there is anything to fill...
  if (!any(is.na(rowSums(pkGrpVal)))) {
    message("No missing peaks present.")
    return(object)
  }

  ## Split the object by file and define the peaks for which
  objectL <- vector("list", length(MSnbase::fileNames(object)))
  pkAreaL <- objectL

  ## We need "only" a list of OnDiskMSnExp, one for each file but
  ## instead of filtering by file we create small objects to keep
  ## memory requirement to a minimum.
  req_fcol  <- MSnbase::requiredFvarLabels("OnDiskMSnExp")
  min_fdata <- object@featureData@data[, req_fcol]
  rt_range  <- range(pkArea[, c("rtmin", "rtmax")])
  if (xcms::hasAdjustedRtime(object))
    min_fdata$retentionTime <- xcms::adjustedRtime(object)
  for (i in 1:length(MSnbase::fileNames(object))) {
    fd <- min_fdata[min_fdata$fileIdx == i, ]
    fd$fileIdx <- 1L
    objectL[[i]] <- new(
      "OnDiskMSnExp",
      processingData = new("MSnProcess",
                           files = MSnbase::fileNames(object)[i]),
      featureData = new("AnnotatedDataFrame", fd),
      phenoData = new("NAnnotatedDataFrame",
                      data.frame(sampleNames = "1")),
      experimentData = new("MIAPE",
                           instrumentManufacturer = "a",
                           instrumentModel = "a",
                           ionSource = "a",
                           analyser = "a",
                           detectorType = "a"))
    ## Extract intensities only for peaks that were not found in a sample.
    pkAreaL[[i]] <- pkArea[is.na(pkGrpVal[, i]), , drop = FALSE]
  }
  rm(pkGrpVal)
  rm(pkArea)
  rm(min_fdata)

  ph <- processHistory(object, type = "Peak detection")
  findPeakMethod <- "unknown"
  mzCenterFun <- "wMean"
  if (length(ph)) {
    if (is(ph[[1]], "XProcessHistory")) {
      prm <- ph[[1]]@param
      findPeakMethod <- .param2string(ph[[1]])
      ## Check if the param class has a mzCenterFun slot
      if (.hasSlot(prm, "mzCenterFun"))
        mzCenterFun <- prm@mzCenterFun
    }
  }
  cp_colnames <- colnames(xcms::chromPeaks(object))

  ## Now rename that to the correct function name in xcms.
  mzCenterFun <- paste("mzCenter", gsub(mzCenterFun, pattern = "mzCenter.", replacement = "", fixed = TRUE), sep=".")
  prepared_data <- list("fdef"=fdef,
                        "mzCenterFun"=mzCenterFun,
                        "object"=object,
                        "pkAreaL"=pkAreaL,
                        "objectL"=objectL,
                        "sampleIndex"=as.list(1:length(objectL)),
                        "cp_colnames"=cp_colnames,
                        "param" = fcp)

  .checkdir(output)
  readr::write_rds(prepared_data, output)
  return()
}


#' Intensive peak extraction step performed file-wise (parallel in cluster workflows)
#'
#' @param input the file index in sample_info
#' @param output the tmp directory for storage of extracted peaks
#' @param index the index of the sample in the sample list (i.e., first file=1, second file=2...)
#' @import xcms
#' @importFrom readr read_rds
#' @importFrom readr write_rds
#' @return a res object for the final peak integration.
#' @export
cFillChromPeaksStep2 <- function(input, output, index){

  # Input: The location of the prepared file from cFillChromPeaksStep1
  # Output: The file-wise output
  prepared_data <- readr::read_rds(input)

  res <-
    .getChromPeakData(
      object      = prepared_data$objectL[[index]],
      peakArea    = prepared_data$pkAreaL[[index]],
      sample_idx  = prepared_data$sampleIndex[[index]],
      cn          = prepared_data$cp_colnames,
      mzCenterFun = prepared_data$mzCenterFun
    )
  # Make separate folder
  .checkdir(output)
  readr::write_rds(res, output)
  return("completed")
}



#' fill chrom peak function 3
#'
#' @param inputFromStep1 Output of cFillChromPeaksStep1() function
#' @param inputFromStep2 Output of cFillChromPeaksStep2() function
#' @param output the output which contains the result of the preprocessing
#'
#' @return a preprocessed xcms object
#' @import xcms
#' @importFrom MSnbase fileNames
#' @importFrom readr read_rds
#' @importFrom readr read_rds
#' @importFrom S4Vectors extractROWS
#' @export
cFillChromPeaksStep3 <- function(inputFromStep1, inputFromStep2, output){

  # Load the data
  prepared_data <- readr::read_rds(inputFromStep1)
  objectL     <- prepared_data$objectL
  pkAreaL     <- prepared_data$pkAreaL
  sample_idx  <- prepared_data$sampleIndex
  cn          <- prepared_data$cp_colnames
  mzCenterFun <- prepared_data$mzCenterFun
  object      <- prepared_data$object
  fdef        <- prepared_data$fdef
  param       <- prepared_data$param
  print(param)
  msLevel     <- 1
  rm(prepared_data)

  res_list <- list()
  for (i in seq_along(inputFromStep2)){
    res_list[[i]] <- readr::read_rds(inputFromStep2[i])
  }

  # START MERGING and save xset_integrated
  res <- do.call(rbind, res_list)
  ## cbind the group_idx column to track the feature/peak group.
  res <- cbind(
    res,
    group_idx = unlist(lapply(pkAreaL, function(z) z[, "group_idx"]), use.names = FALSE)
  )


  ## Remove those without a signal
  res <- res[!is.na(res[, "into"]), , drop = FALSE]
  if (nrow(res) == 0) {
    warning("Could not integrate any signal for the missing ",
            "peaks! Consider increasing 'expandMz' and 'expandRt'.")
    return(object)
  }

  ## Intermediate cleanup of objects.
  rm(pkAreaL)
  gc()

  ## Get the msFeatureData:
  newFd <- new("MsFeatureData")
  newFd@.xData <- .copy_env(object@msFeatureData)
  object@msFeatureData <- new("MsFeatureData")
  incr <- nrow(xcms::chromPeaks(newFd))
  counter = 1
  for (i in unique(res[, "group_idx"])) {
    if (counter %% 100==0){
      print(counter)
    }
    counter = counter+1
    fdef$peakidx[[i]] <- c(fdef$peakidx[[i]],
                           (which(res[, "group_idx"] == i) + incr))
  }
  ## Combine feature data with those from other MS levels
  fdef <- rbind(
    fdef,
    xcms::featureDefinitions(newFd)[xcms::featureDefinitions(newFd)$ms_level != msLevel,,drop = FALSE])

  if (!any(colnames(fdef) == "ms_level")){
    fdef$ms_level <- 1L
  } else {
    fdef <- fdef[order(fdef$ms_level), ]
  }

  ## Define IDs for the new peaks; include fix for issue #347
  maxId <- max(as.numeric(
    sub("M", "", sub("^CP", "", rownames(xcms::chromPeaks(newFd))))))
  if (maxId < 1)
    stop("chromPeaks matrix lacks rownames; please update ",
         "'object' with the 'updateObject' function.")
  toId <- maxId + nrow(res)
  rownames(res) <- sprintf(
    paste0("CP", "%0", ceiling(log10(toId + 1L)), "d"),
    (maxId + 1L):toId)

  chromPeaks(newFd) <- rbind(xcms::chromPeaks(newFd),
                             res[, -ncol(res)])
  cpd           <- S4Vectors::extractROWS(xcms::chromPeakData(newFd), rep(1L, nrow(res)))
  cpd[,]        <- NA
  cpd$ms_level  <- as.integer(msLevel)
  cpd$is_filled <- TRUE
  if (!any(colnames(chromPeakData(newFd)) == "is_filled"))
    xcms::chromPeakData(newFd)$is_filled <- FALSE

  xcms::chromPeakData(newFd) <- rbind(xcms::chromPeakData(newFd), cpd)
  rownames(xcms::chromPeakData(newFd)) <- rownames(xcms::chromPeaks(newFd))
  xcms::featureDefinitions(newFd) <- fdef
  lockEnvironment(newFd, bindings = TRUE)

  object@msFeatureData <- newFd

  ## Add a process history step
  ph <-
    XProcessHistory(
      param = param,
      date. = date(),
      type. = "Missing peak filling",
      fileIndex = 1:length(MSnbase::fileNames(object)),
      msLevel = msLevel
    )

  xset_integrated <- addProcessHistory(object, ph) ## this validates the object.
  readr::write_rds(xset_integrated, output)
  return()
}
