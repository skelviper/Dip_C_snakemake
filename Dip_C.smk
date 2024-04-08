
# Run Dip-C pipeline with snakemake 
# @author zliu
# @date 2022-09-02

#############CONFIG#################

import os

#input
SAMPLES = [i.split(sep='_')[0] for i in os.listdir("./Rawdata")]

configfile: "Dip_C_snakemake/config.yaml"

#############RULE_ALL###############

rule all:
    input:
        expand("processed/{sample}/mapping/{sample}.sort.bam",sample=SAMPLES),
        expand("result/cleaned_pairs/c12/{sample}.pairs.gz",sample=SAMPLES)
        #Hi-C part pairs info
        #expand("result/dip_pairs/{sample}.dip.pairs.gz",sample=SAMPLES),
        #Hi-C part 3d info
        # expand("processed/{sample}/3d_info/{sample}.{res}.align.rms.info",sample=SAMPLES if config["if_structure"] else [],res=["20k","50k","200k","1m"] if config["if_structure"] else []),
        # expand("processed/{sample}/3d_info/{res}.{rep}.3dg", sample=SAMPLES if config["if_structure"] else [],
        #     res=["20k","50k","200k","1m"] if config["if_structure"] else [],
        #     rep=list(range(5)) if config["if_structure"] else []),
        # expand("result/cif_cpg/{sample}.{res}.{rep}.cpg.cif", sample=SAMPLES if config["if_structure"] else [],
        #     res=["20k","50k","200k","1m"] if config["if_structure"] else [],
        #     rep=list(range(5)) if config["if_structure"] else []),
    threads: 20
    shell:"""
         bash ./Dip_C_snakemake/scripts/generateStat.sh
    """


#############END_RULE_ALL##############

# I tried fastqc on a random select file, experiment by Deng YJ, no adapter found.
rule bwa_map:
    input:
        ref_genome=config["refs"][config["ref_genome"]]["bwa_mem2_index"],
        DNA1="Rawdata/{sample}/{sample}_R1.fq.gz",
        DNA2="Rawdata/{sample}/{sample}_R2.fq.gz",
    output:
        bam = "processed/{sample}/mapping/{sample}.sort.bam",
        bamidx = "processed/{sample}/mapping/{sample}.sort.bam.bai"
    threads: config["resources"]["bwa_cpu_threads"],
    resources:
        nodes = config["resources"]["bwa_cpu_threads"],
    params:
        extra=r"-R '@RG\tID:{sample}\tPL:ILLUMINA\tSM:{sample}'",

    shell:"""

        bwa-mem2 mem -5SP -t{threads} {params.extra} {input.ref_genome} {input.DNA1} {input.DNA2} | samtools sort -@{threads} -o {output.bam} -
        samtools index -@{threads} {output.bam} {output.bamidx}

        """

include: "rules/scHiC_2dprocess.rules"
include: "rules/scHiC_3dprocess.rules"
