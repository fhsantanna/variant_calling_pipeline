SAMPLES = ["0288CVL","0548CVL_S2S95","31600CVL","34661CVL","34712CVL","4699CVL","7063CVL","7212CVL","7591CVL"]


rule all:
  input:
    expand("output/variant_call/{sample}.tsv", sample=SAMPLES),
    expand("output/consensus/{sample}_consensus.fa", sample=SAMPLES),
    expand("output/depth/{sample}_primertrim_sorted.depth", sample=SAMPLES)


rule fastp:
    input:
        read_R1="data/reads/{sample}_cat_R1.fastq.gz",
        read_R2="data/reads/{sample}_cat_R2.fastq.gz"
    output:
        fastp_R1="output/fastp/{sample}_R1_fastp.fastq.gz",
        fastp_R2="output/fastp/{sample}_R2_fastp.fastq.gz"
    conda:
        "envs/mapping.yaml"
    shell:
        "fastp -i {input.read_R1} -I {input.read_R2} -o {output.fastp_R1} -O {output.fastp_R2}"

rule index_genome:
    conda:
        "envs/mapping.yaml"
    shell:
        "bwa index data/MN908947.fas"

rule map_reads:
    input:
        "data/MN908947.fas",
        "output/fastp/{sample}_R1_fastp.fastq.gz",
        "output/fastp/{sample}_R2_fastp.fastq.gz"
    output:
        "output/bwa/{sample}_sorted.bam"
    conda:
        "envs/mapping.yaml"    
    shell:
        "bwa mem -t 10 {input} | samtools sort | samtools view -F 4 -o {output}"

rule mark_duplicates:
    input:
        "output/bwa/{sample}_sorted.bam"
    output:
        bam="output/mark_dup/{sample}_markdup.bam",
        metrics="output/mark_dup/{sample}_markdup.txt"
    conda:
        "envs/mapping.yaml"    
    shell:
        "/opt/gatk-4.2.0.0/gatk MarkDuplicates I={input} O={output.bam} M={output.metrics} REMOVE_DUPLICATES=true"

rule trim_primers:
    input:
        "output/mark_dup/{sample}_markdup.bam"
    output:
        "output/primer_trim/{sample}_primertrim.bam"
    conda:
        "envs/mapping.yaml"
    shell:
        "ivar trim -b data/nCoV-2019.primer.bed -p {output} -i {input} -e"

rule sort_bam:
    input:
        "output/primer_trim/{sample}_primertrim.bam"
    output:
        "output/primer_trim/{sample}_primertrim_sorted.bam"
    conda:
        "envs/mapping.yaml"
    shell:
        "samtools sort -o {output} {input}"

rule index_bam:
    input:
        "output/primer_trim/{sample}_primertrim_sorted.bam"
    conda:
        "envs/mapping.yaml"
    shell:
        "samtools index {input}"

rule compute_coverage:
    input:
        "output/primer_trim/{sample}_primertrim_sorted.bam"
    output:
        "output/depth/{sample}_primertrim_sorted.depth"
    conda:
        "envs/mapping.yaml"
    shell:
        "samtools depth -a {input} > {output}"

rule generate_consensus:
    input:
        "output/primer_trim/{sample}_primertrim_sorted.bam"
    output:
        "output/consensus/{sample}_consensus.fa"
    conda:
        "envs/mapping.yaml"
    shell:
        "samtools mpileup -aa -A -d 0 -B -Q 0 {input} | ivar consensus -p {output} -m 5 -n N"

rule variant_call:
    input: 
        "output/primer_trim/{sample}_primertrim_sorted.bam"
    output:
        "output/variant_call/{sample}.tsv"
    conda:
        "envs/mapping.yaml"
    shell:
        "samtools mpileup -A -d 0 -B --reference data/MN908947.fas -Q 0 {input} | ivar variants -p {output} -t 0.03 -m 5 -g data/MN908947.gff -r data/MN908947.fas"

"""
rule variant filtering:
    input: 
        "output/variant_call/{sample}.tsv"
    output:
        "output/variant_call/{sample}_filtered.tsv"
    conda:
        "envs/mapping.yaml"
    shell:
        "ivar filtervariants -p {output} {input}"
"""
