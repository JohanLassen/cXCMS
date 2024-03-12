from gwf import Workflow, AnonymousTarget
import os
import re

gwf = Workflow()


### Helper functions 
def runR(rfile, inputfile, outputfile, cores = 1, memory="12g", walltime="01:00:00"):
    inputs = [inputfile]
    outputs = [outputfile]
    options = {
        "nodes": 1,
        'cores': cores,
        'memory': memory,
        "walltime":walltime,
        "account" : "forensics"
    }
    spec = '''
    echo jobinfo $SLURM_JOBID
    Rscript {rfile} {inputfile} {outputfile}
    '''.format(rfile = rfile, inputfile=inputfile, outputfile=outputfile)
    return AnonymousTarget(inputs=inputs, outputs=outputs, spec=spec, options=options)

def runR_integration(rfile, inputfile, outputfile, extra, cores = 1, memory="8g", walltime="01:00:00"):
    inputs = [inputfile]
    outputs = [outputfile]
    options = {
        "nodes": 1,
        'cores': cores,
        'memory': memory,
        "walltime":walltime,
        "account" : "forensics"
    }
    spec = '''
    echo jobinfo $SLURM_JOBID
    Rscript {rfile} {inputfile} {outputfile} {extra}
    '''.format(rfile = rfile, inputfile=inputfile, outputfile=outputfile, extra = extra)
    return AnonymousTarget(inputs=inputs, outputs=outputs, spec=spec, options=options)



def filter_strings(strings, substring):
    return [s for s in strings if substring in s]


def list_files_recursive(path):
    """Lists all files in a directory and its subdirectories recursively."""
    files = []
    for root, _, filenames in os.walk(path):
        for filename in filenames:
            file_path = os.path.join(root, filename)
            files.append(file_path)
    return files


### Identify files
all_files = list_files_recursive("/faststorage/project/forensics/Kirstine/Kathrine_GHB_mzml")
mzml_files = filter_strings(all_files, "mzML")
#pos_files = filter_strings(mzml_files, "Data_to_XCMS_pos")

os.makedirs("./tmp/", exist_ok=True)


### Generate directory structure
tmp_folder = "./tmp/"
peaks_folder = tmp_folder+"peak_picked/"
integration_folder = tmp_folder+"peak_integrated/"
results_folder = "./results/"

os.makedirs(tmp_folder, exist_ok=True)
os.makedirs(peaks_folder, exist_ok=True)
os.makedirs(integration_folder, exist_ok=True)
os.makedirs(results_folder, exist_ok=True)

#################### Step 1: Peak Picking

inputs = mzml_files
outputs = [re.sub(".*Kirstine/", "", filnavn) for filnavn in inputs]
outputs = [re.sub("/", "_", filnavn) for filnavn in outputs]
outputs = [peaks_folder+re.sub("mzML", "rds", filnavn) for filnavn in outputs]
#print(outputs[:10])
# Iterate through all files
index = 0
for inputfile, outputfile in zip(inputs, outputs):
    gwf.target_from_template(
        "peak_picking"+str(index),
        runR_integration(
            rfile="./scripts/step1_peakcalling.R",
            inputfile=inputfile, 
            outputfile=outputfile,
            extra="./sample_overview.csv")
    )
    index += 1
    #break


#################### Step 2: Merging files

output_merged = "./tmp/peak_picked_merge_dataset.rds"
inputs_path = peaks_folder

gwf.target_from_template(
    "peak_merging",
    runR(rfile="./scripts/step2_merging.R",
        inputfile=inputs_path, 
        outputfile=output_merged)
)

#################### Step 3: Grouping

output_group1 = "./tmp/peak_grouped1_dataset.rds"
input = output_merged

gwf.target_from_template(
    "peak_grouping1",
    runR(rfile="./scripts/step3_grouping1.R",
        inputfile=input, 
        outputfile=output_group1)
)

#################### Step 4: Alignment

output_alignment= "./tmp/peak_aligned_dataset.rds"
input = output_group1

gwf.target_from_template(
    "peak_alignment",
    runR(rfile="./scripts/step4_alignment.R",
        inputfile=input, 
        outputfile=output_alignment)
)

#################### Step 5: Grouping2

output_group2 = "./tmp/peak_grouped2_dataset.rds"
input = output_alignment

gwf.target_from_template(
    "peak_grouping2",
    runR(rfile="./scripts/step5_grouping2.R",
        inputfile=input, 
        outputfile=output_group2)
)

#################### Step 6: Integration1

output_integration1 = "./tmp/peak_integration1_dataset.rds"
input = output_group2

gwf.target_from_template(
    "peak_integration1",
    runR(rfile="./scripts/step6_integration1.R",
        inputfile=input, 
        outputfile=output_integration1)
)

#################### Step 7: Peak integration 2

output_peakint2 = [re.sub("peak_picked", "peak_integrated", x) for x in outputs]
input = output_integration1
print(input)

index = 1 # R uses indexing from 1.
for outputfile in output_peakint2:
    gwf.target_from_template(
        "peak_integration2_"+str(index),
        runR_integration(rfile="./scripts/step7_integration2.R",
            inputfile=input, 
            outputfile=outputfile, 
            memory = "56g",
            extra = index)
    )
    index += 1
    #break


#################### Step 8: Integration3

output_integration3 = "./tmp/peak_integration3_dataset.rds"
input = integration_folder

gwf.target_from_template(
    "peak_integration3",
    runR_integration(rfile="./scripts/step8_integration3.R",
        inputfile=input, 
        outputfile=output_integration3,
        extra = output_integration1
        )
)

#################### Step 9: Generate names for features in output
output_feature_names = "./results/output_feature_names_pos.csv"
input = output_integration3

gwf.target_from_template(
    "output_feature_names_pos",
    runR(rfile="./scripts/step9_writeoutput.R",
        inputfile=input, 
        outputfile=output_feature_names
        )
)

#################### Step 10: Generate CAMERA output
output_feature_names = "./results/diffrep.csv"
input = output_integration3

gwf.target_from_template(
    "isotope_annotation",
    runR(rfile="./scripts/step10_adduct_annotation.R",
        inputfile=input, 
        outputfile=output_feature_names,
        memory="32g",
        walltime="10:00:00"
        )
)

