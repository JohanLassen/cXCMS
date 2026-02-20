



#' Create output dir
#'
#' This function is a utility function that creates the output directory if it does not already exist
#'
#' @param output the output path
#'
#' @return None
.checkdir <- function(output){
  dir_path <- dirname(output)
  if (!file.exists(dir_path)){
    dir.create(dir_path, recursive=TRUE)
  }
  return("dir created")
}




.mse_combine <- function(x) {

  res <- x[[1L]]
  res@sampleData <- do.call(MsCoreUtils::rbindFill, lapply(x, MsExperiment::sampleData))
  res@spectra <- do.call(c, lapply(x, ProtGenerics::spectra))
  sl <- lapply(x, function(z) z@sampleDataLinks[["spectra"]])
  nsamp <- lengths(x)
  nsamp <- c(0, cumsum(nsamp)[-length(nsamp)])
  nspec <- vapply(sl, nrow, NA_integer_)
  nspec <- c(0, cumsum(nspec)[-length(nspec)])
  res@sampleDataLinks[["spectra"]] <- do.call(
    rbind, mapply(function(z, i, j) {
      z[, 1L] <- z[, 1L] + i
      z[, 2L] <- z[, 2L] + j
      z
    }, sl, nsamp, nspec, SIMPLIFY = FALSE, USE.NAMES = FALSE))
  res
}


.xmse_combine <-
  function(x) {
    res <- .mse_combine(x)
    nsamp <- lengths(x)
    nsamp <- c(0, cumsum(nsamp)[-length(nsamp)])
    res@chromPeaks <- do.call(MsCoreUtils::rbindFill, mapply(function(z, i) {
      z <- xcms::chromPeaks(z)
      z[, "sample"] <- z[, "sample"] + i
      z
    }, x, nsamp, SIMPLIFY = FALSE, USE.NAMES = FALSE))
    rownames(res@chromPeaks) <- xcms:::.featureIDs(nrow(res@chromPeaks), "CP")
    res@chromPeakData <- do.call(
      MsCoreUtils::rbindFill,
      lapply(x, xcms::chromPeakData, return.type = "data.frame"))
    rownames(res@chromPeakData) <- rownames(res@chromPeaks)
    res@processHistory <- do.call(c, lapply(x, xcms::processHistory))
    res
  }


#' Single file peak identification for parallel processing and/or RAM optimization
#'
#' @param input filename of the mzML input files
#' @param output output filename
#' @param cwp native xcms CentWaveParam() function output
#' @param metadata sample metadata
#' @param settings settings list from settings.yaml
#' @importFrom xcms findChromPeaks MergeNeighboringPeaksParam refineChromPeaks filterRt
#' @importFrom MsExperiment readMsExperiment sampleData
#' @importFrom Spectra MsBackendMzR
#' @importFrom ProtGenerics filterSpectra spectra
#' @importFrom MsCoreUtils rbindFill
#' @importFrom readr write_rds
#' @importFrom tibble tibble
#' @return peak picked files in the save_folder location (./tmp/peaks_identified)
#' @export
cFindChromPeaksStep1 <- function(input, output, cwp, metadata = NULL, settings = NULL){

  if (typeof(input)!="character"){stop("Input must be a string\n")}
  if (!grepl("mzML|mzXML", input)){stop("Input must be a mzML or mzXML file\n")}
  if (is.null(metadata)){
    warning("No groups provided (Experimental groups). Using random group assignment\n")
    metadata <- tibble(
      sample_name = gsub("[.]mzML", "", basename(input)),
      group = sample(paste0("group", 1:2), 1)
    )}
  if (is.null(metadata$sample_group)){
    warning("No groups provided (Experimental groups). Using random group assignment\n")
    metadata$groups = sample(paste0("group", 1:2), 1)
  }

  # Get RT filter from settings or use defaults
  rt_min <- if (!is.null(settings$general$rt_filter_min)) settings$general$rt_filter_min else 30

rt_max <- if (!is.null(settings$general$rt_filter_max)) settings$general$rt_filter_max else 600

  object <-
    MsExperiment::readMsExperiment(
      spectraFiles = input,
      backend = Spectra::MsBackendMzR(),
      sampleData = metadata
    )
  object <- ProtGenerics::filterSpectra(object, filterRt, rt = c(rt_min, rt_max))

  object  <- xcms::findChromPeaks(object, cwp, msLevel = 1L)

  mnpp  <- MergeNeighboringPeaksParam(expandRt = 5)
  object <- refineChromPeaks(object, param = mnpp, chunkSize = 1)

  readr::write_rds(object, output)
}



#' Fast concatenation of XCMS peak picked files
#'
#' @param inputs the peak picked files
#' @param output the msnbase element of peak picked files
#' @import xcms
#' @importFrom readr read_rds write_rds
#' @importFrom purrr map
#' @return an XcmsExperiment ready for peak grouping and alignment
#' @export
cFindChromPeaksStep2 <- function(inputs, output){
  ## This function concatenates the files by building a string followed by concatenation.
  ## By doing it this way we reduce running time from n^2 to n allowing us to concatenate several thousand files in minutes rather than days.

  if (length(inputs)<2){
    stop("Less than two files in the inputs. Please use *all* peak called files as a vector in the inputs parameter")
  }

  object <- inputs |> purrr::map(~readr::read_rds(.x), .progress=TRUE) |> .xmse_combine()
  readr::write_rds(object, output)
  return(object)
}





