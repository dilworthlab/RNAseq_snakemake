import pandas as pd
import numpy as np
import glob
import os
import yaml
from pathlib import Path


from snakemake.io import expand, glob_wildcards
from snakemake.utils import min_version
from snakemake.logging import logger


# Minimum snakemake version
min_version("5.24.1")

# configuration file
configfile:"./config/config.yaml"

# Tabular configuration
samples = pd.read_csv(config["sample_info"], "\t").set_index("samples")


# Multiqc configuration
home = str(Path.home())
datadict = {'log_filesize_limit': 2000000000}

files = glob.glob(os.path.join(home, ".*.yaml"))
if os.path.join(home, ".multiqc_config.yaml") in files:
    print('~/.multiqc_config.yaml is present')
    with open(os.path.join(home, ".multiqc_config.yaml"), 'r') as file:
        values = yaml.safe_load(file)

        if 'log_filesize_limit' in str(values):
            print('log_filesize_limit is already set')
        if values is None:
            with open(os.path.join(home, ".multiqc_config.yaml"), 'w') as file:
                docs = yaml.dump(datadict, file)
                print('added log_filesize_limit to multiqc log file')
        else:
            with open(os.path.join(home, ".multiqc_config.yaml"), 'r') as file:
                new_yaml = yaml.safe_load(file)
                print(type(new_yaml))
                new_yaml.update(datadict)


            with open(os.path.join(home, ".multiqc_config.yaml"),'w') as file:
                yaml.safe_dump(new_yaml, file)
                print('log_filesize_limit: 2000000000 is set')

else:
    with open(os.path.join(home, ".multiqc_config.yaml"), 'w') as file:
        documents = yaml.dump(datadict, file)
        print('made new .multiqc_config.yaml')






# Getting directory containing rawreads
wdir = os.getcwd()
for filepath, dirs, allfiles in os.walk(wdir):
    for file_ in allfiles:
        if file_.endswith(".fastq.gz") or file_.endswith(".fastq"):
            READS_DIR = filepath
            ROOT_DIR = os.path.relpath(filepath, wdir)

# fastq or fasta.gz
readfiles = os.listdir(READS_DIR)
filename = readfiles[1]
SUFFIX = filename.split(".")[1]


logger.info(f'This is the RawReads dir: {READS_DIR}')

# Describing wildcards
SINGLE_READ = f'{READS_DIR}/{{fastqfile}}_{{read}}.{SUFFIX}'


READS = set(glob_wildcards(SINGLE_READ).read)
FASTQFILES = set(glob_wildcards(SINGLE_READ).fastqfile)

logger.info(f'These are the read wildcards: {READS}')
logger.info(f'These are the fastqfile wildcards: {FASTQFILES}')



rule all:
    input:
        'logs/Salmon/Mus_musculus_index/make_index.log',
        expand("Analysis_Results/QC_Rawreads/{fastqfile}_{read}_fastqc.html", fastqfile=FASTQFILES, read=READS),
        expand("Analysis_Results/QC_Rawreads/{fastqfile}_{read}_fastqc.zip", fastqfile=FASTQFILES, read=READS),
        "Analysis_Results/QC_Rawreads/Rawreads_QC.html"




#Quality Control FastQC
rule QCrawreads_Fastqc:
    input:
        f'{ROOT_DIR}/{{fastqfile}}_{{read}}.{SUFFIX}',
    output:
        ('Analysis_Results/QC_Rawreads/{fastqfile}_{read}_fastqc.html'),
        ('Analysis_Results/QC_Rawreads/{fastqfile}_{read}_fastqc.zip')
    log:
        'logs/fastqc_rawreads/{fastqfile}_{read}.log'
    threads: 6
    resources:
        mem_mb=6000,
        runtime=1440
    shell:
        """
        fastqc {input[0]} --outdir=./Analysis_Results/QC_Rawreads &>> {log}
        """

rule Compileresults_QC:
    input: expand('Analysis_Results/QC_Rawreads/{fastqfile}_{read}_fastqc.html', read=READS, fastqfile=FASTQFILES)
    output:
     html="Analysis_Results/QC_Rawreads/Rawreads_QC.html"
    log:
        'logs/compileresults/QC.log'
    shell:
        "multiqc ./Analysis_Results/QC_Rawreads --force -v -o ./Analysis_Results/QC_Rawreads -n Rawreads_QC.html &>> {log} "



#Organism = config['Organism']
#if Organism_index

rule GetIndex_Salmon:
    input:
        transcriptome='Mus_musculus/gencode.vM25.transcripts.fa.gz'
    output:
        'logs/Salmon/Mus_musculus_index/make_index.log'
    resources:
        mem_mb=40000,
        runtime=2160
    shell:
        "salmon index -t {input.transcriptome} -i Mus_musculus/salmon_index &>> {output} "



rule Quantify_Salmon:
    input:
        'logs/Salmon/Mus_musculus_index/make_index.log',
        read1=os.path.join(ROOT_DIR, ('{fastqfile}_1.' + SUFFIX)),
        read2=os.path.join(ROOT_DIR, ('{fastqfile}_2.' + SUFFIX))
    output:
        quants=directory('quants/{fastqfile}_quant')
    log:
        'logs/Salmon/Quantify/{fastqfile}.log'
    params:
        index='Mus_musculus/salmon_index'
    threads: 8
    resources:
        mem_mb=40000,
        runtime=2160
    shell:
        " salmon quant -i Mus_musculus/salmon_index -l A -1 {input.read1} -2 {input.read2} -p {threads} --validateMappings -o {output.quants} &>> {log} "
