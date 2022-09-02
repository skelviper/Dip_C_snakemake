
# Run Dip-C pipeline with snakemake 
# @author zliu
# @date 2022-09-02

#############CONFIG#################

import os

#input
SAMPLES = [i.split(sep='_')[0] for i in os.listdir("./Rawdata")]
#SPLIT = ["atac","ct"]
#SAMPLES = os.listdir("./Rawdata")
#SAMPLES = ["E752001"]

configfile: "Dip_C_snakemake/config.yaml"

#############RULE_ALL###############

rule all:
    input:
        #Hi-C part pairs info
        expand("result/cleaned_pairs/c12/{sample}.pairs.gz",sample=SAMPLES),
        expand("result/dip_pairs/{sample}.dip.pairs.gz",sample=SAMPLES),
        #Hi-C part 3d info
        expand("processed/{sample}/3d_info/{sample}.{res}.align.rms.info",sample=SAMPLES if config["if_structure"] else [],res=["20k","50k","200k","1m"] if config["if_structure"] else []),
        expand("processed/{sample}/3d_info/{res}.{rep}.3dg", sample=SAMPLES if config["if_structure"] else [],
            res=["20k","50k","200k","1m"] if config["if_structure"] else [],
            rep=list(range(5)) if config["if_structure"] else []),
        expand("result/cif_cpg/{sample}.{res}.{rep}.cpg.cif", sample=SAMPLES if config["if_structure"] else [],
            res=["20k","50k","200k","1m"] if config["if_structure"] else [],
            rep=list(range(5)) if config["if_structure"] else []),


#############END_RULE_ALL##############

# I tried fastqc on a random select file, experiment by Deng YJ, no adapter found.
# rule trim_adapter:
#     input:
#         R1 = "Rawdata/{smaple}/{sample}_R1.fastq.gz",
#         R2 = "Rawdata/{smaple}/{sample}_R2.fastq.gz",
#     output:
#         R1 = "processed/{sample}/DNA/{sample}_R1.fastq.gz",
#         R2 = "processed/{sample}/DNA/{sample}_R2.fastq.gz",
#     shell:"""
#         source ~/miniconda3/etc/profile.d/conda.sh
#         conda activate py3
#         set -u

#         bwa-mem2 mem -5SP -t{threads} {params.extra} {input.ref_genome} {input.DNA1} {input.DNA2} | samtools sort -@{threads} -o {output.bam} -
#         samtools index -@{threads} {output.bam} {output.bamidx}

#         set +u
#         conda deactivate
#     """

include: "rules/scHiC_2dprocess.rules"
include: "rules/scHiC_3dprocess.rules"
