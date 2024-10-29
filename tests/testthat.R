# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html

#library(testthat)
#library(cXCMS)

# test_check("cXCMS")
# files <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
# out1  <- c("test1.rds", "test2.rds")
#
# # Write tests using the README examples
#
# test_that("cFindChromPeaksStep1", {
#   # Test that the function runs without errors
#
#   expect_warning(cFindChromPeaksStep1(files[1], "test.rds", cwp = xcms::CentWaveParam()))
# })
#
# #
# # test_that("cFindChromPeaksStep2", {
# #   # Test that the function runs without errors
# #   for (i in seq_along(files)){
# #     cFindChromPeaksStep1(files[i], out1[i], cwp = xcms::CentWaveParam())
# #   }
# #   expect_message(cFindChromPeaksStep2(out1, "test.rds"))
# # })
