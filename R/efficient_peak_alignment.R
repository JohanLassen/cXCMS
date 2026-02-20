#' @importFrom dplyr group_split
#' @importFrom purrr imap map
.applyRtAdjToChromPeaks <- function(x, rtraw, rtadj) {
  if (!is.list(rtraw) | !is.list(rtadj))
    stop("'rtraw' and 'rtadj' are supposed to be lists!")
  if (length(rtraw) != length(rtadj))
    stop("'rtraw' and 'rtadj' have to have the same length!")

  cols_to_adjust <- c("rt", "rtmin", "rtmax")
  samples <- dplyr::group_split(data.frame(x[,cols_to_adjust]), x[,"sample"], .keep = FALSE)

  adj <- samples |>
    purrr::imap(~{
      xcms:::.applyRtAdjustment(
        x = as.matrix(.x),
        rtraw = rtraw[[.y]],
        rtadj = rtadj[[.y]]
      )},
      .progress = TRUE) |>
    purrr::map(matrix, ncol = length(cols_to_adjust)) |>
    purrr::map(~{colnames(.x) <- cols_to_adjust; .x})
  x[,cols_to_adjust] <- do.call(rbind, adj)
  x
}


setMethod(
  "adjustRtime", signature(object = "MsExperiment",
                           param = "PeakGroupsParam"),
  function(object, param, msLevel = 1L, ...) {
    if (!nrow(xcms:::peakGroupsMatrix(param))) {
      xcms:::peakGroupsMatrix(param) <- xcms:::adjustRtimePeakGroups(object, param)
    }
    fidx <- as.factor(fromFile(object))
    rt_raw <- split(rtime(object), fidx)
    rt_adj <- xcms:::.adjustRtime_peakGroupsMatrix(
      rt_raw, xcms:::peakGroupsMatrix(param), smooth = param@smooth,
      span = param@span, family = param@family,
      subset = param@subset, subsetAdjust = param@subsetAdjust)
    pt <- vapply(object@processHistory, xcms:::processType, character(1))
    idx_pg <- xcms:::.match_last(xcms:::.PROCSTEP.PEAK.GROUPING, pt, nomatch = -1L)
    if (idx_pg > 0) {ph <- object@processHistory[idx_pg]} else {ph <- list()}
    object <- xcms::dropFeatureDefinitions(object)
    object@spectra@backend@spectraData$rtime_adjusted <- unlist(rt_adj, use.names = FALSE)
    object@chromPeaks <- .applyRtAdjToChromPeaks(
      xcms:::.chromPeaks(object), rtraw = rt_raw, rtadj = rt_adj)
    xph <- xcms:::XProcessHistory(
      param = param, type. = xcms:::.PROCSTEP.RTIME.CORRECTION,
      fileIndex = seq_along(object), msLevel = msLevel)
    phist <- object@processHistory
    phist <- c(phist, ph, list(xph))
    object@processHistory <- phist
    object
  })
  