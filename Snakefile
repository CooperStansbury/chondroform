from datetime import datetime
import pandas as pd
import yaml
from pathlib import Path
import re
import os
import sys
from tabulate import tabulate


BASE_DIR = Path(workflow.basedir)
configfile: str(BASE_DIR) + "/config/config.yaml"

# big picture variables
OUTPUT = config['output_path']
print("\nOUTPUT PATH:")
print(OUTPUT)

# load input paths
input_path = os.path.abspath(config['inputs'])
input_df = pd.read_csv(input_path, comment="#")
samples = input_df['sample_id'].to_list()

# get new path names
output_paths = [OUTPUT + "raw_fastq/" + x + '.fastq.gz' for x in samples]

print("OUTPUT PATHS:")
print(output_paths)

# print statements
print("\n----- CONFIG VALUES -----")
for key, value in config.items():
    print(f"{key}: {value}")
    
    
print("\n----- INPUT VALUES -----")
print(
    tabulate(
        input_df, 
        headers='keys', 
        tablefmt='psql',
        showindex=False,
    )
)


rule all:
    input:
        expand(OUTPUT + "raw_fastq/{sid}.fastq.gz", sid=samples),
        OUTPUT + 'references/annotations.gtf',
        OUTPUT + 'references/reference.fa',
        OUTPUT + 'references/reference.mmi',
        expand(OUTPUT + "mapping/{sid}.bam.bai", sid=samples),
        expand(OUTPUT + "count_matrix/{sid}.counts.tsv", sid=samples),
     

rule get_gtf:
    input:
        config['gtf_path']
    output:
        OUTPUT + 'references/annotations.gtf.gz',
    shell:
        """cp {input} {output} """
        
        
rule prep_annotations:
    input:
        refgenome=OUTPUT + 'references/annotations.gtf.gz',
    output:
        ref=OUTPUT + 'references/annotations.gtf',
        flag=touch(OUTPUT + 'references/annotations.done')
    shell:
        "cat {input} | gzip -d > {output.ref}"
        
        
rule get_reference:
    input:
        config['ref_path']
    output:
        OUTPUT + 'references/reference.fa.gz',
    shell:
        """cp {input} {output} """
        
    
rule prep_reference:
    input:
        refgenome=OUTPUT + 'references/reference.fa.gz'
    output:
        ref=OUTPUT + 'references/reference.fa',
        flag=touch(OUTPUT + 'references/reference.done')
    shell:
        "cat {input} | gzip -d > {output.ref}"
        
        
rule minimap2_index:
    input:
        refgenome=OUTPUT + 'references/reference.fa.gz'
    output:
        OUTPUT + 'references/reference.mmi'
    threads:
        config['threads']
    conda:
        "aligner"
    shell:
        "minimap2 -t {threads} -d {output} {input.refgenome}"
        
        
rule get_fastq:
    input:
        input_df['file_path'].to_list()
    output:
        output_paths
    run:
        from shutil import copyfile
        for i, refPath in enumerate(input):

            outPath = output[i]
            copyfile(refPath, outPath)
            

rule align_reads:
    input:
        fastq=OUTPUT + "raw_fastq/{sid}.fastq.gz",
        ref=OUTPUT + 'references/reference.fa.gz',
        refindex=OUTPUT + 'references/reference.mmi',
    output:        
        bam=OUTPUT + 'mapping/{sid}.bam',
    params:
        args=config['minimap2_args'],
    threads:
        int(config['threads'])
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    log:
        OUTPUT + "mapping/{sid}.log",
    conda:
        "aligner"
    shell:
        """minimap2 {params.args} -t {threads} \
        {input.ref} {input.fastq} | samtools sort \
        -@ {threads} -O bam -o {output.bam} """
        
        
        
rule samtools_index:
    input:
        OUTPUT + 'mapping/{sid}.bam'
    output:
        OUTPUT + 'mapping/{sid}.bam.bai'
    conda:
        "aligner"
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    threads:
        int(config['threads'])
    conda:
        "aligner"
    shell:
        """samtools index -@ {threads} {input}"""



rule htseq_count:
    input:
        bam=OUTPUT + 'mapping/{sid}.bam',
        annotations=OUTPUT + 'references/annotations.gtf',
    output:
        OUTPUT + "count_matrix/{sid}.counts.tsv"
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    conda:
        "count_matrix"
    params:
        FEATURE_TYPE='gene',
        IDATTR='gene_id',
        ADDITIONAL_ATTRIBUTES='gene_name',
    threads:
        config['threads']
    shell:
        """htseq-count {input.bam}  {input.annotations} \
        -t {params.FEATURE_TYPE} -i {params.IDATTR} \
        --additional-attr {params.ADDITIONAL_ATTRIBUTES} \
        --with-header -n  {threads} > {output}"""
       



# rule build_db:
#     input:
#         OUTPUT + 'references/annotations.gtf',
#     output:
#         OUTPUT + "isoquant/annotations.db",
#     conda:
#         "../envs/isoquant.yml"
#     shell:
#         """python scripts/build_isoquant_db.py {input} {output}"""
# 
# 
# rule run_isoquant:
#     input:
#         ref=OUTPUT + 'references/reference.fa',
#         db=OUTPUT + "isoquant/annotations.db",
#         bam=OUTPUT + "mapping/{sid}.tagged.bam",
#         bam_index=OUTPUT + "mapping/{sid}.tagged.bam.bai",
#     output:
#         directory(OUTPUT + "isoquant/{sid}"),
#         touch(OUTPUT + "isoquant/{sid}.done"),
#     params:
#         outdir=OUTPUT + "isoquant",
#         prefix=lambda wildcards: wildcards.sid
#     conda:
#         "../envs/isoquant.yml"
#     threads:
#         config['threads']
#     wildcard_constraints:
#         sid='|'.join([re.escape(x) for x in set(samples)]),
#     shell:
#         """isoquant.py --reference {input.ref} --genedb {input.db} \
#         --threads {threads} --count_exons --prefix {params.prefix} \
#         --bam {input.bam} --data_type 'nanopore' \
#         --bam_tags 'CB,UB,RD' \
#         --complete_genedb -o {params.outdir}"""