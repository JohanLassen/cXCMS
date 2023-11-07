



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
  return()
}


#' Single file peak identification for parallel processing and/or RAM optimization
#'
#' @param input filename of the mzML input files
#' @param output output filename
#' @param cwp native xcms CentWaveParam() function output
#' @param groups the group variable used in new("NAnnotatedDataFrame", data.frame("sample_name"=sample_name, "group"=groups)). The group variable might bias the analysis and in big datasets it might beneficial setting groups = None to do random group assignment
#' @importFrom xcms findChromPeaks
#' @importFrom MSnbase readMSData
#' @return peak picked files in the save_folder location (./tmp/peaks_identified)
#' @export
cFindChromPeaksStep1 <- function(input, output, cwp, groups = NULL){

  if (typeof(input)!="character"){stop("Input must be a string\n")}
  if (!grepl("mzML|mzXML", input)){stop("Input must be a mzML or mzXML file\n")}
  if (!grepl("[.]Rdata|[.]rdata", output)){stop("Output file must be a rdata file type (.Rdata or .rdata). Just edit the name of the output files\n")}
  if (is.null(groups)){
    warning("No groups provided (Experimental groups). Using random group assignment\n")
    groups = sample(paste0("group", 1:2), 1)}

  sample_name <- gsub("[.].*|.*[/]", "", input)
  # Return if output file already has been preprocessed
  if (file.exists(output)) return()

  # Make the XCMSnExp object
  loaded_file <-
    MSnbase::readMSData(
      files = input,
      pdata = new("NAnnotatedDataFrame", data.frame("sample_name"=sample_name, "group"=groups)),
      mode = "onDisk",
      msLevel. = 1L
    )

  # Call peaks
  peaks_file  <- xcms::findChromPeaks(loaded_file, cwp)
  # Save result
  save(peaks_file, file=output)
}



#' Fast concatenation of XCMS peak picked files
#'
#' @param inputs the peak picked files
#' @param output the msnbase element of peak picked files
#' @import xcms
#' @importFrom readr read_rds
#' @importFrom readr write_rds
#' @return an XCMSset ready for peak grouping and alignment
#' @export
cFindChromPeaksStep2 <- function(inputs, output){
  ## This function concatenates the files by building a string followed by concatenation.
  ## By doing it this way we reduce running time from n^2 to n allowing us to concatenate several thousand files in minutes rather than days.

  if (length(inputs)<2){
    stop("Less than two files in the inputs. Please use *all* peak called files as a vector in the inputs parameter")}

  # First element of string
  concatenate_string <- "c("
  for (i in 1:length(inputs)){

    # arbitrary name for the concatenated list
    name <- paste0("peak", i)
    file <- inputs[i]
    peaks_file <- readr::read_rds(file) # loads XCMSnExp object w. identified peaks
    assign(name, peaks_file)

    # Prepare for another file...
    if (i<length(inputs)){concatenate_string <- paste0(c(concatenate_string, name, ","), collapse = "")}
    # Or end the concatenation function c()
    if (i==length(inputs)){concatenate_string <- paste0(c(concatenate_string,  name, ")"), collapse = "")}
  }

  # Evaluate string
  object <- eval(parse(text = concatenate_string))
  readr::write_rds(object, output)
}







