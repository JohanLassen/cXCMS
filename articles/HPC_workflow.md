# HPC workflow

The goal of cXCMS is to provide the framework for efficient high
performance computing (HPC) untargeted LCMS data pre-processing. Also,
the package memory optimizes the workflow in laptops environments (MacOS
or Windows), allowing preprocessing of larger sample sizes on small PCs.
The source code has been used for the project [Large-Scale metabolomics:
Predicting biological age using 10,133 routine untargeted LC-MS
measurements](NA).

## HCP Introduction

This tutorial is meant for high performance computing users on linux
systems. The tutorial walks through the following parts:

- Prerequisites
- Installation
- Workflow Setup
- Run workflow
- Details

#### Prerequisites

Before starting make sure you have
[anaconda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)
and [git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)
installed.

Currently, the workflow is based on [gwf](https://gwf.app/) and has been
tested on Slurm although SGE and LSF also should work. We invite
snakemake users to push a workflow as well.

Finally, it is important you have the data available as mzML files.

#### Installation

Open the terminal and navigate to your project directory:

    cd ~/home/users/path/to/experiment/

Download the workflow from github

    git clone http://github.com/JohanLassen/xcms_workflow

Navigate to the workflow folder

    cd xcms_workflow

Now you should install the dependencies using conda. This might take
some time.

    conda env create -n xcms_workflow -f environment.yml

Activate the conda environment

    conda activate xcms_workflow

#### Workflow Setup

Before running the workflow, edit the setup.yaml file. The setup file
has two categories of parameters:

1.  Those making the code work

- run_id (needed): Identifier of your run. Smart if both positive and
  negative ionization must be run. Will be used to name submitted jobs.

- input_path (needed): The absolute path to the input data.
  Subdirectories are allowed and will be included in the file identifier
  in the output. E.g., the data directory could be ./data/positive and
  contain the subdirectories ./data/positive/batch1 and
  ./data/positive/batch2

- output_path (strongly recommended): The output directory. If it does
  not exist the workflow will generate it automatically.

- sample_overview_path (optional): Currently on experimental stage. Used
  for defining groups of the samples. The file must be semicolon
  separated with a column named “sample_group” and a column named
  “sample_name”. The sample_names must be (at least) unique substrings
  of the filenames, e.g., for file ./datafolder/sample_id_123.mzML the
  sample_name “sample_id_123” is sufficient.

- CAMERA_RULES (optional): path to the adduct annotation rules (CAMERA R
  package). Default ones are included in the workflow.

2.  Those improving your xcms results:

- As indicated the workflow does peak picking → grouping → alignment →
  grouping and finally integration of missing peaks. The parameters of
  each step must be set accordingly to the experiment.

Natively it looks like this and is meant for positive ionization data,
see parameter descriptions below:

    general:
      run_id: run1
      input_path: ./data/
      output_path: ./experiment_results/
      sample_overview_path: ./overview_pos.csv
      CAMERA_RULES: ./camera_rules/rules_pos.csv
    xcms_parameters:
      centwave:
        mzdiff: 0.01
        peakwidth: c(4,25)
        ppm: 12
        prefilter: c(3,500)
        snthresh: 6
      peak_grouping1:
        binSize: 0.01
        bw: 2.5
        maxFeatures: 50
        minFraction: 0.3
      alignment:
        extraPeaks: 5
        family: gaussian
        minFraction: 0.9
        smooth: loess
        span: 0.6
      peak_grouping2:
        binSize: 0.01
        bw: 2.5
        maxFeatures: 50
        minFraction: 0.6

#### Run the workflow

To run the workflow one must initiate gwf (remember to activate the
conda environment). This only has to be done once per workflow
directory.

    gwf config set backend slurm

Now you are ready to submit the jobs. Be sure that the setup file is
correct before - the workflow will submit a bunch of jobs (2\*number of
files).

    gwf run

#### Details

The workflow *generates a schedule file* (.csv file) at the output
directory. The workflow will only do this once per run_id, meaning that
if any of the files are corrupt these can be removed from the schedule
file and the workflow can be resubmitted without recomputing the already
finished parts. If the schedule file is deleted the next workflow
resubmission generates a new one including all original files.

The workflow saves a *copy of the setup file*, meaning you always can
recover the settings that the experiment was run under. Just copy paste
the file to the workflow directory to rerun the analysis.

The tmp folder can be deleted after the run, but it might be beneficial
to keep the peak_picked files if more batches are added later on.
Additionally, it is possible to load the temporary files into an
interactive R session if any preprocessed step needs some quality
assurance, plotting, or another output is desired by the user.

The *outputs* include a feature table with one column describing the
file names (including sub directories separated by underscores) and the
remaining columns representing features. It also returns a peak table
containing adduct annotations. Most xcms users will be familiar with the
formats.
