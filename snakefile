from __future__ import print_function
import os
import fnmatch

SAMPLES = os.listdir("data/libraries")

# PROGGENY_SPL = fnmatch.filter(SAMPLES, 'F*')
# PARENT_SPL = fnmatch.filter(SAMPLES, 'S*')

rule all:
    input:
        expand("data/libraries/{sample}/{sample}.sam", sample=SAMPLES),
        expand("data/libraries/{sample}/{sample}.bam", sample=SAMPLES),
        expand("data/libraries/{sample}/{sample}_sorted.bam", sample=SAMPLES),
        expand("data/libraries/{sample}/{sample}_sorted.bam.bai", sample=SAMPLES),
        expand("data/libraries/{sample}/{sample}_sorted_MD.bam", sample=SAMPLES),
        expand("data/libraries/{sample}/{sample}_sorted_MD.log", sample=SAMPLES),
        expand("data/libraries/{sample}/{sample}_sorted_MD.bam.bai", sample=SAMPLES),
        expand("data/libraries/{sample}/{sample}_sorted_MD.grp", sample=SAMPLES),
        expand("data/libraries/{sample}/{sample}_sorted_MD_recal.bam", sample=SAMPLES),
        expand("data/libraries/{sample}/{sample}_sorted_MD_recal.bam.bai", sample=SAMPLES),
        expand("data/libraries/{sample}/{sample}_sorted_MD_recal.flagstat", sample=SAMPLES)

rule alignment:
    input:
        read1="data/libraries/{sample}/{sample}_R1.fastq.gz",
        read2="data/libraries/{sample}/{sample}_R2.fastq.gz",
        genome="data/genome/schistosoma_mansoni.PRJEA36577.WBPS14.genomic.fa",
    output:
        "data/libraries/{sample}/{sample}.sam"
    params:
        rg=r"@RG\tID:{sample}\tPL:illumina\tLB:{sample}\tSM:{sample}"
    shell:
        'sleep $[ ( $RANDOM % 100 )  + 1 ]s ; '
        'bwa mem -t $(nproc) -M -R "{params.rg}" "{input.genome}" "{input.read1}" "{input.read2}" > "{output}"'

rule bam:
    input:
        "data/libraries/{sample}/{sample}.sam"
    output:
        "data/libraries/{sample}/{sample}.bam"
    shell:
        'sleep 1m $[ ( $RANDOM % 200 )  + 1 ]s ; '
        'samtools view -b -h -@ $(nproc) "{input}" > "{output}"'

rule sorting:
    input:
        "data/libraries/{sample}/{sample}.bam"
    output:
        "data/libraries/{sample}/{sample}_sorted.bam"
    shell:
        'sleep 15s $[ ( $RANDOM % 100 )  + 1 ]s ; '
        'samtools sort -@ $(nproc) -f "{input}" "{output}"'

rule indexing1:
    input:
        "data/libraries/{sample}/{sample}_sorted.bam"
    output:
        "data/libraries/{sample}/{sample}_sorted.bam.bai",
    shell:
        'sleep 1m $[ ( $RANDOM % 200 )  + 1 ]s ; '
        'samtools index "{input}"'

rule mark_duplicates:
    input:
        "data/libraries/{sample}/{sample}_sorted.bam"
    output:
        bam="data/libraries/{sample}/{sample}_sorted_MD.bam",
        log="data/libraries/{sample}/{sample}_sorted_MD.log"
    shell:
        'gatk --java-options "-Xmx2g" MarkDuplicates -I "{input}" -O "{output.bam}" -M "{output.log}" --VALIDATION_STRINGENCY LENIENT --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP $(ulimit -n)'

rule indexing2:
    input:
        "data/libraries/{sample}/{sample}_sorted_MD.bam"
    output:
        "data/libraries/{sample}/{sample}_sorted_MD.bam.bai"
    shell:
        'sleep 1m $[ ( $RANDOM % 200 )  + 1 ]s ; '
        'samtools index "{input}"'

rule bsqr:
    input:
        bam="data/libraries/{sample}/{sample}_sorted_MD.bam",
        genome="data/genome/schistosoma_mansoni.PRJEA36577.WBPS14.genomic.fa",
        sites="data/genome/sm_dbSNP_v7.vcf"
    output:
        table="data/libraries/{sample}/{sample}_sorted_MD.grp",
        bam="data/libraries/{sample}/{sample}_sorted_MD_recal.bam"
    params:
        table=r"data/libraries/{sample}/{sample}_sorted_MD.grp"
    shell:
        'gatk --java-options "-Xmx2g" BaseRecalibrator -R "{input.genome}" -I "{input.bam}" --known-sites "{input.sites}" -O "{output.table}" &&\
         gatk --java-options "-Xmx2g" ApplyBQSR -R "{input.genome}" -I "{input.bam}" --bqsr-recal-file "{params.table}" -O "{output.bam}"'

rule indexing3:
    input:
        "data/libraries/{sample}/{sample}_sorted_MD_recal.bam"
    output:
        "data/libraries/{sample}/{sample}_sorted_MD_recal.bam.bai"
    shell:
        'sleep 1m $[ ( $RANDOM % 200 )  + 1 ]s ; '
        'samtools index "{input}"'

rule stats:
    input:
        "data/libraries/{sample}/{sample}_sorted_MD_recal.bam"
    output:
        "data/libraries/{sample}/{sample}_sorted_MD_recal.flagstat"
    params:
        spl=r"data/libraries/{sample}/{sample}_sorted_MD_recal.bam"
    shell:
        'echo -e "File: {params.spl}" > "{output}" ; \
         samtools flagstat "{input}" >> "{output}" ; \
         echo  "\n" >> "{output}"'

# rule clenup:
#     input:
#         sam="data/libraries/{sample}/{sample}.sam"
#         bam1="data/libraries/{sample}/{sample}_sorted.bam"
#         bam2="data/libraries/{sample}/{sample}_sorted_MD.bam"
#         table="data/libraries/{sample}/{sample}_sorted_MD.grp"
#     shell:
#         'rm "{input.sam}"* "{input.bam1}"* "{input.bam2}"* "{input.table}"'

