from __future__ import print_function
import os
import fnmatch

SAMPLES = os.listdir("data/libraries")

# PROGGENY_SPL = fnmatch.filter(SAMPLES, 'F*')
# PARENT_SPL = fnmatch.filter(SAMPLES, 'S*')

rule all:
    input:
        expand("data/libraries/{sample}/{sample}_sorted.bam", sample=SAMPLES),
        expand("data/libraries/{sample}/{sample}_sorted.bam.bai", sample=SAMPLES),
        expand("data/libraries/{sample}/{sample}_sorted_MD.bam", sample=SAMPLES),
        expand("data/libraries/{sample}/{sample}_sorted_MD.log", sample=SAMPLES),
        expand("data/libraries/{sample}/{sample}_sorted_MD.bam.bai", sample=SAMPLES),
        expand("data/libraries/{sample}/{sample}_sorted_MD.grp", sample=SAMPLES),
        expand("data/libraries/{sample}/{sample}_sorted_MD_recal.bam", sample=SAMPLES),
        expand("data/libraries/{sample}/{sample}_sorted_MD_recal.bam.bai", sample=SAMPLES),
        expand("data/libraries/{sample}/{sample}_sorted_MD_recal.flagstat", sample=SAMPLES),
        expand("data/libraries/{sample}/{sample}.gvcf.gz", sample=SAMPLES),
        "data/calling/cerc_prod.gvcf.gz",
        "data/calling/cerc_prod.vcf.gz"

rule alignment:
    input:
        read1="data/libraries/{sample}/{sample}_R1.fastq.gz",
        read2="data/libraries/{sample}/{sample}_R2.fastq.gz",
        genome="data/genome/schistosoma_mansoni.PRJEA36577.WBPS14.genomic.fa"
    output:
        temp("data/libraries/{sample}/{sample}_sorted.bam")
    params:
        rg=r"@RG\tID:{sample}\tPL:illumina\tLB:{sample}\tSM:{sample}"
    shell:
        'sleep $[ ( $RANDOM % 100 )  + 1 ]s ; '
        'bwa mem -t $(nproc) -M -R "{params.rg}" "{input.genome}" "{input.read1}" "{input.read2}" | samtools sort -@8 -o "{output}" -'

rule indexing1:
    input:
        "data/libraries/{sample}/{sample}_sorted.bam"
    output:
        temp("data/libraries/{sample}/{sample}_sorted.bam.bai")
    shell:
        'sleep 1m $[ ( $RANDOM % 200 )  + 1 ]s ; '
        'samtools index "{input}"'

rule mark_duplicates:
    input:
        bam="data/libraries/{sample}/{sample}_sorted.bam",
        bai="data/libraries/{sample}/{sample}_sorted.bam.bai"
    output:
        bam=temp("data/libraries/{sample}/{sample}_sorted_MD.bam"),
        log=protected("data/libraries/{sample}/{sample}_sorted_MD.log")
    shell:
        'gatk --java-options "-Xmx2g" MarkDuplicates -I "{input.bam}" -O "{output.bam}" -M "{output.log}" --VALIDATION_STRINGENCY LENIENT --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP $(ulimit -n)'

rule indexing2:
    input:
        "data/libraries/{sample}/{sample}_sorted_MD.bam"
    output:
        temp("data/libraries/{sample}/{sample}_sorted_MD.bam.bai")
    shell:
        'sleep 1m $[ ( $RANDOM % 200 )  + 1 ]s ; '
        'samtools index "{input}"'

rule bsqr:
    input:
        bam="data/libraries/{sample}/{sample}_sorted_MD.bam",
        bai="data/libraries/{sample}/{sample}_sorted_MD.bam.bai",
        genome="data/genome/schistosoma_mansoni.PRJEA36577.WBPS14.genomic.fa",
        sites="data/genome/sm_dbSNP_v7.vcf"
    output:
        table=temp("data/libraries/{sample}/{sample}_sorted_MD.grp"),
        bam=protected("data/libraries/{sample}/{sample}_sorted_MD_recal.bam")
    params:
        table=r"data/libraries/{sample}/{sample}_sorted_MD.grp"
    shell:
        'gatk --java-options "-Xmx2g" BaseRecalibrator -R "{input.genome}" -I "{input.bam}" --known-sites "{input.sites}" -O "{output.table}" &&\
         gatk --java-options "-Xmx2g" ApplyBQSR -R "{input.genome}" -I "{input.bam}" --bqsr-recal-file "{params.table}" -O "{output.bam}"'

rule indexing3:
    input:
        "data/libraries/{sample}/{sample}_sorted_MD_recal.bam"
    output:
        protected("data/libraries/{sample}/{sample}_sorted_MD_recal.bam.bai")
    shell:
        'sleep 1m $[ ( $RANDOM % 200 )  + 1 ]s ; '
        'samtools index "{input}"'

rule stats:
    input:
        bam="data/libraries/{sample}/{sample}_sorted_MD_recal.bam",
        bai="data/libraries/{sample}/{sample}_sorted_MD_recal.bam.bai"
    output:
        protected("data/libraries/{sample}/{sample}_sorted_MD_recal.flagstat")
    params:
        spl=r"data/libraries/{sample}/{sample}_sorted_MD_recal.bam"
    shell:
        'echo -e "File: {params.spl}" > "{output}" ; \
         samtools flagstat "{input.bam}" >> "{output}" ; \
         echo  "\n" >> "{output}"'

rule calling:
    input:
        bam="data/libraries/{sample}/{sample}_sorted_MD_recal.bam",
        bai="data/libraries/{sample}/{sample}_sorted_MD_recal.bam.bai",
        genome="data/genome/schistosoma_mansoni.PRJEA36577.WBPS14.genomic.fa",
        sites="data/genome/sm_dbSNP_v7.vcf"
    output:
        "data/libraries/{sample}/{sample}.gvcf.gz"
    shell:
        'sleep 30s $[ ( $RANDOM % 100 )  + 1 ]s ; '
        'gatk --java-options "-Xmx2g" HaplotypeCaller -R "{input.genome}" -I "{input.bam}" -D "{input.sites}" --output-mode EMIT_ALL_ACTIVE_SITES -ERC GVCF -O "{output}"'

rule combining:
    input:
        gvcfs=expand("data/libraries/{sample}/{sample}.gvcf.gz", sample=SAMPLES),
        genome="data/genome/schistosoma_mansoni.PRJEA36577.WBPS14.genomic.fa",
        sites="data/genome/sm_dbSNP_v7.vcf"
    output:
        temp("data/calling/cerc_prod.gvcf.gz")
    run:
        gvcfs=" --variant ".join(input.gvcfs)
        shell('sleep 30s $[ ( $RANDOM % 200 )  + 1 ]s ; '
        'gatk --java-options "-Xmx9g" CombineGVCFs -R "{input.genome}" --variant {gvcfs} -D "{input.sites}" -O "{output}"')

rule merging:
    input:
        gvcf="data/calling/cerc_prod.gvcf.gz",
        genome="data/genome/schistosoma_mansoni.PRJEA36577.WBPS14.genomic.fa",
        sites="data/genome/sm_dbSNP_v7.vcf"
    output:
        protected("data/calling/cerc_prod.vcf.gz")
    shell:
        'sleep 30s $[ ( $RANDOM % 200 )  + 1 ]s ; '
        'gatk --java-options "-Xmx2g"  GenotypeGVCFs -R "{input.genome}" -V "{input.gvcf}" -D "{input.sites}" -O "{output}"'
