# Single file peak identification for parallel processing and/or RAM optimization

Single file peak identification for parallel processing and/or RAM
optimization

## Usage

``` r
cFindChromPeaksStep1(input, output, cwp, metadata = NULL, settings = NULL)
```

## Arguments

- input:

  filename of the mzML input files

- output:

  output filename

- cwp:

  native xcms CentWaveParam() function output

- metadata:

  sample metadata

- settings:

  settings list from settings.yaml

## Value

peak picked files in the save_folder location (./tmp/peaks_identified)
