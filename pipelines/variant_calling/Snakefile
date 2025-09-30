# Snakefile 
# Prerequisites (on Gadi): source /g/data/pq84/malaria/snakemake_pipeline/snakemake_setup.sh


################################################ Define paths using config file #########################################################################


configfile: "config/config.yaml"
fasta_path=config['fasta_path']
picard_path=config['picard_path']
gatk_path=config['gatk_path']
known_indels_path=config['known_indels_path']
known_sites_path=config['known_sites_path']
bed_file_path=config['bed_file_path']
source_dir=config['source_dir']


################################################ Set up: define input files, wildcards, etc. #########################################################################


import pandas as pd
import math
import re

## Defining the samples to be used for the {sample} wildcard
SAMPLES, = glob_wildcards(f'{source_dir}/{{sample}}_1.fastq.gz') 

## Chromosome names for bed files and subsetting jobs
with open(config['bed_file_path']) as f:
    CHROMOSOME = f.read().splitlines()
    CHROMOSOME = [p.split('\t')[0] for p in CHROMOSOME]
## Create bed for each chromosome/contig - to be used for sub-setting jobs (eg joint calling)

with open(config['bed_file_path']) as f:
    chrom_bed = f.read().splitlines()
    chrom_bed = [p.split('\t') for p in chrom_bed]
    
chrom_bed_df = pd.DataFrame(chrom_bed, columns = ("chrom", "start", "end"))
chrom_bed_df = chrom_bed_df.astype({'start': 'int32', 'end': 'int32'})
chrom_bed_df_small = chrom_bed_df[chrom_bed_df['end'] < 100000]
chrom_bed_df = chrom_bed_df[chrom_bed_df['end'] > 100000]

## End columns
for i in range(1, 11):
    col_name = 'seg_end' + str(i)
    chrom_bed_df[col_name] = chrom_bed_df.apply(lambda x: math.floor(x['end']/10 * i), axis=1, result_type = 'expand')

## Start columns
columns = [col for col in chrom_bed_df if 'seg_end' in col]
for col in columns:
    col_str = re.sub("end.*", "start", col)
    col_num = col.replace("seg_end", "")
    col_num = int(col_num) + 1
    col_name = col_str + str(col_num)
    chrom_bed_df[col_name] = chrom_bed_df.apply(lambda x: x[col] + 1, axis = 1, result_type = 'expand')

# Wrangle into appropriate format and shift start position to 1
chrom_bed_df['start'] = chrom_bed_df.apply(lambda x: x['start'] + 1, axis = 1)
chrom_bed_df_small['start'] = chrom_bed_df_small.apply(lambda x: x['start'] + 1, axis = 1)

chrom_bed_df = chrom_bed_df.rename(columns={'start': 'seg_start1'}).drop(['end', 'seg_start11'], axis = 1)
chrom_bed_df = pd.wide_to_long(chrom_bed_df, stubnames = ['seg_end', 'seg_start'], i = 'chrom', j = 'segment').reset_index(level = ['chrom']) # seg_end tells the fucntion to use the end of this as the new row values for j (segment)

chrom_bed_df_small = chrom_bed_df_small.rename(columns = {'start': 'seg_start', 'end': 'seg_end'})
chrom_bed_df = pd.concat([chrom_bed_df, chrom_bed_df_small])

# Create numpy arrays
CHROMOSOME_INTERVALS = chrom_bed_df['chrom'] + ':' + chrom_bed_df['seg_start'].astype(str) + '-' + chrom_bed_df['seg_end'].astype(str)
CHROMOSOME_INTERVALS = CHROMOSOME_INTERVALS.to_numpy()


################################################## Define final files #######################################################################


rule all:
    input:
        "output/calling/consensus/Consensus.vcf.gz"


###################################################### Rules ##############################################################################


# Define local rules - not run with scheduler
localrules: all, bam_input_list

# Concatenate chromosome-based consensus VCFs 
rule concat_vcfs:
    input:
        expand("output/calling/consensus/{chromosome}_consensus.vcf.gz", chromosome = CHROMOSOME)
    params:
        vcf = lambda w: " ".join(expand("output/calling/consensus/{chromosome}_consensus.vcf.gz", chromosome = CHROMOSOME))
    output:
        vcf="output/calling/consensus/Consensus.vcf.gz",
        tbi="output/calling/consensus/Consensus.vcf.gz.tbi"
    shell:
        """
        bcftools concat -o {output.vcf} {input}
        bcftools index -t -o {output.tbi} {output.vcf}
        """

# Take a consenus of GATK and bcftools
rule consensus_of_vcfs:
    input:
        bcftools="output/calling/bcftools/bcftools_genotyped_{chromosome}.vcf.gz",
        gatk="output/calling/gatk/joint/gatk_genotyped_{chromosome}.vcf.gz"
    output:
        txt=temp("output/calling/consensus/{chromosome}.txt"),
        vcf=temp("output/calling/consensus/{chromosome}_consensus.vcf.gz"),
        tbi=temp("output/calling/consensus/{chromosome}_consensus.vcf.gz.tbi")
    params:
        header = lambda w: f"'%CHROM\\t%POS\\n'"
    shell:
        """
        bcftools query -f {params.header} {input.bcftools} > {output.txt}
        bcftools filter -R {output.txt} -o {output.vcf} {input.gatk}
        bcftools index -t -o {output.tbi} {output.vcf}
        """

# Concatenate chromosome-based consensus VCFs 
def get_intervals_by_chromosome(wildcards):
    return [f"output/calling/bcftools/intervals/bcftools_genotyped_intervals_{intervals}.vcf.gz"
            for intervals in CHROMOSOME_INTERVALS if intervals.startswith(wildcards.chromosome+':')]

rule concat_bcftools:
    input:
        get_intervals_by_chromosome
    output:
        vcf=temp("output/calling/bcftools/bcftools_genotyped_{chromosome}.vcf.gz"),
        tbi=temp("output/calling/bcftools/bcftools_genotyped_{chromosome}.vcf.gz.tbi")
    params:
        interval_list=lambda w, input: " ".join(input)
    shell:
        """
        bcftools concat -o {output.vcf} {params.interval_list}
        bcftools index -t -o {output.tbi} {output.vcf}
        """

# Run bcftools for each chromosome interval
rule bcftools_caller:
    input:
        bam=expand("output/bam_recal/{sample}_recalibrated.bam", sample = SAMPLES),
        input_bam_files="output/calling/bcftools/input_bam_files.list",
        fasta=fasta_path
    output:
        vcf="output/calling/bcftools/intervals/bcftools_genotyped_intervals_{intervals}.vcf.gz",
        tbi="output/calling/bcftools/intervals/bcftools_genotyped_intervals_{intervals}.vcf.gz.tbi"
    params:
        bed=lambda w: w.intervals
    shell:
        """
        bcftools mpileup --threads 2 -f {input.fasta} -b {input.input_bam_files} -r {params.bed} | bcftools call --threads 2 -m -Oz -a FORMAT/GQ,FORMAT/GP,INFO/PV4 -v -o {output.vcf}
        bcftools index --threads 2 -t -o {output.tbi} {output.vcf}
        """
    
# Create input list of bam files for bcftools

rule bam_input_list:
    input:
        expand("output/bam_recal/{sample}_recalibrated.bam", sample = SAMPLES)
    output:
        temp("output/calling/bcftools/input_bam_files.list")
    run:
        import glob

        bam_list = glob.glob('output/bam_recal/*_recalibrated.bam')

        file = open('output/calling/bcftools/input_bam_files.list', 'w')
        for item in bam_list:
            file.write(item+"\n")

        file.close()

# Joint-call variants

rule joint_genotyping:
    input:
        vcf="output/calling/gatk/gvcf/GATK_combined.g.vcf.gz",
        fasta=fasta_path
    output:
        "output/calling/gatk/joint/gatk_genotyped_{chromosome}.vcf.gz"
    params:
        gatk=gatk_path,
        bed=lambda w: w.chromosome
    shell:
        """
        java -Djava.iodir=1000m -Xms3200m -Xmx3600m -jar {params.gatk} \
        -T GenotypeGVCFs \
        -nt 3 \
        -R {input.fasta} \
        -L {params.bed} \
        -V {input.vcf} \
        -o {output}
        """

# Combine gVCFs

rule combine_gvcfs:
    input:
        fasta=fasta_path,
        vcf=expand("output/calling/gatk/gvcf/{sample}.g.vcf.gz", sample = SAMPLES)
    output:
        temp("output/calling/gatk/gvcf/GATK_combined.g.vcf.gz")
    params:
        gatk=gatk_path,
        gvcfs=lambda w: " -V " + " -V ".join(expand("output/calling/gatk/gvcf/{sample}.g.vcf.gz", sample = SAMPLES))
    shell:
        """
        java -Djava.iodir=1000m -Xms3200m -Xmx3600m -jar {params.gatk} \
        -T CombineGVCFs \
        -R {input.fasta} \
        {params.gvcfs} \
        -o {output}
        """

# Call haplotypes - GATK

rule haplotype_caller:
    input:
        bam="output/bam_recal/{sample}_recalibrated.bam",
        fasta=fasta_path
    output:
        "output/calling/gatk/gvcf/{sample}.g.vcf.gz"
    params:
        gatk=gatk_path
    shell:
        """
        java -Djava.iodir=1000m -Xms3200m -Xmx3600m -jar {params.gatk} \
        -T HaplotypeCaller \
        -ERC GVCF \
        --minPruning 3 \
        --maxNumHaplotypesInPopulation 200 \
        --max_alternate_alleles 3 \
        --variant_index_type LINEAR \
        --variant_index_parameter 128000 \
        -contamination 0.0 \
        -G Standard \
        -R {input.fasta} \
        -I {input.bam} \
        -o {output}
        """

# Get recalibrated bams

rule print_reads:
    input:
        bam="output/bam_recal/{sample}_dupmarked_realigned.bam",
        table="output/bam_recal/{sample}_dupmarked_realigned_recal.table",
        fasta=fasta_path,
        bed=bed_file_path
    output:
        "output/bam_recal/{sample}_recalibrated.bam"
    params:
        gatk=gatk_path
    shell:
        """
        java -Djava.iodir=1000m -Xms3200m -Xmx3600m -jar {params.gatk} \
        -T PrintReads \
        -R {input.fasta} \
        --intervals {input.bed} \
        -I {input.bam} \
        -BQSR {input.table} \
        -o {output}
        """

# Create table for recalibration

rule table_for_base_recal:
    input:
        bam="output/bam_recal/{sample}_dupmarked_realigned.bam",
        fasta=fasta_path,
        sites=known_sites_path,
        bed=bed_file_path
    output:
        temp("output/bam_recal/{sample}_dupmarked_realigned_recal.table")
    params:
        gatk=gatk_path
    shell:
        """
        java -Djava.iodir=1000m -Xms3200m -Xmx3600m -jar {params.gatk} \
        -T BaseRecalibrator \
        -R {input.fasta} \
        -I {input.bam} \
        --intervals {input.bed} \
        -knownSites {input.sites} \
        -o {output}
        """

# Realign BAM files by indels

rule indel_realigner:
    input:
        bam="output/bam_recal/{sample}_dupmarked_reheader.bam",
        fasta=fasta_path,
        indels=known_indels_path,
        bed=bed_file_path,
        targets="output/bam_recal/{sample}_dupmarked_realigner.intervals"
    output:
       temp("output/bam_recal/{sample}_dupmarked_realigned.bam")
    params:
        gatk=gatk_path
    shell:
        """ 
        java -Djava.iodir=1000m -Xms3200m -Xmx3600m -jar {params.gatk} \
        -T IndelRealigner \
        --consensusDeterminationModel KNOWNS_ONLY \
        -LOD 0.4 \
        -R {input.fasta} \
        -I {input.bam} \
        --intervals {input.bed} \
        -known {input.indels} \
        -targetIntervals {input.targets} \
        -o {output}
        """

# Create Realignment Targets

rule realigner_target_creator:
    threads: 5
    input:
        bam="output/bam_recal/{sample}_dupmarked_reheader.bam",
        fasta=fasta_path,
        indels=known_indels_path,
        bed=bed_file_path
    output:
        temp("output/bam_recal/{sample}_dupmarked_realigner.intervals")
    params:
        gatk=gatk_path
    shell:
        """
        java -Djava.iodir=1000m -Xms3200m -Xmx3600m -jar {params.gatk} \
        -T RealignerTargetCreator \
        -nt {threads} \
        -R {input.fasta} \
        -I {input.bam} \
        --intervals {input.bed} \
        -known {input.indels} \
        -o {output}
        """

# Update headers and index bam files

rule update_header_and_index:
    input:
        "output/bam_recal/{sample}_dupmarked.bam"
    output:
        bam_output=temp("output/bam_recal/{sample}_dupmarked_reheader.bam"),
        bam_index=temp("output/bam_recal/{sample}_dupmarked_reheader.bam.bai")
    params:
        header=lambda w: f"'s,^@RG.*,@RG\\tID:{w.sample}\\tSM:{w.sample}\\tLB:None\\tPL:Illumina,g'"
    shell:
        """
        samtools view -H {input} | \
        sed {params.header} | \
        samtools reheader - {input} > {output.bam_output}

        samtools index {output.bam_output} 
        """

# Mark duplicates

rule mark_duplicates:
    input: 
        bam="output/sorted/{sample}_sorted.bam",
        bam_index="output/sorted/{sample}_sorted.bam.bai"
    output: 
        dup_marked=temp("output/bam_recal/{sample}_dupmarked.bam"),
        metrics_file=temp("output/bam_recal/{sample}_picard_metrics_file.txt")
    params:
        picard=picard_path
    shell:
        """
        java -Djava.iodir=1000m -Xms3200m -Xmx3600m -jar {params.picard} \
        MarkDuplicates AS=TRUE VALIDATION_STRINGENCY=LENIENT \
        I={input.bam} \
        O={output.dup_marked} \
        M={output.metrics_file}
        """
        
# Sort BAM for recalibration & index

rule index_sorted_bam: 
    input: 
        "output/sorted/{sample}_sorted.bam"
    output:
        temp("output/sorted/{sample}_sorted.bam.bai") 
    shell:
        "samtools index {input}"

rule samtools_sort:
    threads: 5
    input: 
        "output/mapped_reads/{sample}.bam"
    output:
        temp("output/sorted/{sample}_sorted.bam")
    shell:
        "samtools sort -@ {threads} {input} > {output}"

# Map reads 

rule bwa_map:
    threads: 5
    input:
        fasta_path,
        f"{source_dir}/{{sample}}_1.fastq.gz",
        f"{source_dir}/{{sample}}_2.fastq.gz"
    output:
        "output/mapped_reads/{sample}.bam"
    params:
        header=lambda w: f"@RG\\\\tID:{w.sample}\\\\tPL:ILLUMINA"
    shell:
         "bwa mem -t {threads} -M -R {params.header} {input} | samtools view -u -S - | samtools sort -n -o {output}"
