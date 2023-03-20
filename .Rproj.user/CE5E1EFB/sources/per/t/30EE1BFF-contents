


#' Single file peak identification for parallel processing and/or RAM optimization
#'
#' @param sample_info data.frame (or tibble) containing the columns: "filename" and "group"
#' @param cwp native xcms CentWaveParam() function output
#' @param index which file should be run in xcms
#' @param save_folder the folder for saving the peak picked files. It is recommended only to change this if running custom workflows.
#' @import xcms
#' @return peak picked files in the save_folder location (./tmp/peaks_identified)
#' @export
cFindChromPeaks <- function(sample_info, cwp, index = 1, save_folder = NULL){

  if (!is.data.frame(sample_info)){
    stop("'sample_info' should be a dataframe containing absolute mzml file locations and group names!")}
  if (!("filename" %in% colnames(sample_info))){
    stop("'sample_info' should contain the column 'filename' containing absolute mzml file locations!\nExample: /home/user/Documents/lcms_project/data/fileA.mzML")}
  if (!("group" %in% colnames(sample_info))){
    warning("'sample_info' should contain the column 'group' containing the groups of the experiment. None found, using random group assignment")}

  # generate sample_names consistent throughout the full workflow
  sample_info$sample_name <-  gsub(".*/", "", sample_info$filename)

  # Make sure directory exists
  if (is.null(save_folder)){
    save_folder = "./tmp/peaks_identified/"
    if (!dir.exists("./tmp")){
      dir.create("./tmp")
      dir.create("./tmp/peaks_identified")
    }
  }

  # Return if output file already has been preprocessed
  output_file <- paste0(save_folder, gsub("[.].*", ".Rdata", sample_info$sample_name[index]))
  print(output_file)
  if (file.exists(output_file)) return()

  # Make the XCMSnExp object
  loaded_file <-
    MSnbase::readMSData(
      files = sample_info$filename[index],
      pdata = new("NAnnotatedDataFrame", sample_info[, c("sample_name", "group")]),
      mode = "onDisk",
      msLevel. = 1L
    )
  # Call peaks
  peaks_file  <- xcms::findChromPeaks(loaded_file, cwp)
  # Save result
  save(peaks_file, file=output_file)
}



#' Fast concatenation of XCMS peak picked files
#'
#' @param peaks_identified_path the location of the peak picked files
#' @import xcms
#' @return an XCMSset ready for peak grouping and alignment
#' @export
collect_files <- function(peaks_identified_path = "./tmp/peaks_identified/"){
  ## This function concatenates the files by building a string followed by concatenation.
  ## By doing it this way we reduce running time from n^2 to n allowing us to concatenate several thousand files in minutes rather than days.

  files <- list.files(peaks_identified_path, full.names = T)
  if (length(files)==0){
    stop(paste0("No peak called files in ", peaks_identified_path, ". The path is wrong or files haven't been peak picked yet"))}

  # First element of string
  concatenate_string <- "c("
  for (i in 1:length(files)){
    # arbitrary name for the concatenated list
    name <- paste0("peak", i)
    file <- files[i]
    load(file) # loads XCMSnExp object w. identified peaks
    assign(name, peaks_file)

    # Prepare for another file...
    if (i<length(files)){concatenate_string <- paste0(c(concatenate_string, name, ","), collapse = "")}
    # Or end the concatenation function c()
    if (i==length(files)){concatenate_string <- paste0(c(concatenate_string,  name, ")"), collapse = "")}
  }

  # Evaluate string
  xset_peaksIdentified <- eval(parse(text = concatenate_string))
  save(xset_peaksIdentified, file = "./tmp/xset_peaksIdentified.Rdata")
  return(xset_peaksIdentified)
}







