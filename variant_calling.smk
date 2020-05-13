from __future__ import print_function
import os
import fnmatch

SAMPLES = os.listdir("data/libraries")

rule all:
    input:
        expand("data/libraries/{sample}/{sample}.gvcf.gz", sample=SAMPLES),
        "data/calling/cerc_prod.gvcf.gz",
        "data/calling/cerc_prod.vcf.gz"

rule calling:
    input:
        bam="data/libraries/{sample}/{sample}_sorted_MD_recal.bam",
        genome="data/genome/schistosoma_mansoni.PRJEA36577.WBPS14.genomic.fa",
        sites="data/genome/sm_dbSNP_v7.vcf"
    output:
        "data/libraries/{sample}/{sample}.gvcf.gz"
    shell:
        'sleep 30s $[ ( $RANDOM % 100 )  + 1 ]s ; '
        'gatk --java-options "-Xmx2g" HaplotypeCaller -R "{input.genome}" -I "{input.bam}" -D "{input.sites}" --output-mode EMIT_ALL_ACTIVE_SITES -ERC GVCF -O "{output}"'

rule combining:
    input:
        # gvcfs="data/libraries/{sample}/{sample}.gvcf.gz",
        gvcfs=expand("data/libraries/{sample}/{sample}.gvcf.gz", sample=SAMPLES),
        genome="data/genome/schistosoma_mansoni.PRJEA36577.WBPS14.genomic.fa",
        sites="data/genome/sm_dbSNP_v7.vcf"
    output:
        "data/calling/cerc_prod.gvcf.gz"
    run:
        gvcfs=" --variant ".join(input.gvcfs)
        shell('sleep 30s $[ ( $RANDOM % 200 )  + 1 ]s ; '
        'gatk --java-options "-Xmx2g" CombineGVCFs -R "{input.genome}" {gvcfs} -D "{input.sites}" -O "{output}"')

rule merging:
    input:
        gvcf="data/calling/cerc_prod.gvcf.gz",
        genome="data/genome/schistosoma_mansoni.PRJEA36577.WBPS14.genomic.fa",
        sites="data/genome/sm_dbSNP_v7.vcf"
    output:
        "data/calling/cerc_prod.vcf.gz"
    shell:
        'sleep 30s $[ ( $RANDOM % 200 )  + 1 ]s ; '
        'gatk --java-options "-Xmx2g"  GenotypeGVCFs -R "{input.genome}" -V "{input.gvcf}" -D "{input.sites}" -O "{output}"'
