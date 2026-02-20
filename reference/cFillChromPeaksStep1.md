# Preparation for the computational intensive peak filling step

cPrepFillChromPeaks, cGetChromPeakData, and cFillChromPeaks are small
constructs of fillChromPeaks.

By running the functions as suggested in the vignette, we minimize the
RAM usage, thus making the workflow scalable to tens of thousands of
files.

## Usage

``` r
cFillChromPeaksStep1(input, output, param, interval, settings)
```

## Arguments

- input:

  path to XcmsExperiment as returned by the peak alignment step

- output:

  path(s) of the output file(s). This file goes into
  cFillChromPeaksStep2()

- param:

  the object returned by xcms::FillChromPeaksParam()

- interval:

  numeric vector of length 2 with start and end sample indices

- settings:

  settings list from settings.yaml

## Value

a list object with all the necessary objects for the cGetChromPeakData
and cFillChromPeaks functions
