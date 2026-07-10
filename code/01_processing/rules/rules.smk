include: 'common.smk'


rule fastqc:
    input:
        os.path.join(config['dir']['raw'], "{sample}_{num}.fq.gz")
    output:
        html = os.path.join(config['dir']['qc'], "{sample}_{num}_fastqc.html"),
        zip = os.path.join(config['dir']['qc'], "{sample}_{num}_fastqc.zip") # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: "--quiet"
    log:
        os.path.join(config['dir']['log'], "fastqc_{sample}_{num}.log")
    threads: 1
    wrapper:
        "0.79.0/bio/fastqc"


rule cutadapt_5end:
    input:
        os.path.join(config['dir']['raw'], '{sample}_1.fq.gz'),
        os.path.join(config['dir']['raw'], '{sample}_2.fq.gz')
    output:
        fastq1 = os.path.join(config['dir']['results']['fq'], '{sample}_5end_1.fq.gz'),
        fastq2 = os.path.join(config['dir']['results']['fq'], '{sample}_5end_2.fq.gz'),
        qc = os.path.join(config['dir']['qc'], "{sample}_cutadapt_5end.txt")
    params:
        adapters = "-g file:{R1} -G file:{R2}".format(R1=config['adapters']['R1']['5end'],
                                                      R2=config['adapters']['R2']['5end']),
        extra = "--discard-untrimmed",
    log:
        os.path.join(config['dir']['log'], "{sample}_cutadapt_5end.log")
    threads:
        config['threads']
    wrapper:
        "0.79.0/bio/cutadapt/pe"


rule cutadapt_5end_rename_R1:
    input:
        os.path.join(config['dir']['raw'], '{sample}_1.fq.gz'),
        os.path.join(config['dir']['raw'], '{sample}_2.fq.gz')
    output:
        fastq1 = os.path.join(config['dir']['results']['fq'], '{sample}_5end_1.fq.gz'),
        fastq2 = os.path.join(config['dir']['results']['fq'], '{sample}_5end_2.fq.gz'),
        qc = os.path.join(config['dir']['qc'], "{sample}_cutadapt_5end.txt")
    params:
        adapters = "-g file:{R1} -G file:{R2}".format(R1=config['adapters']['R1']['5end'],
                                                      R2=config['adapters']['R2']['5end']),
        extra = lambda wildcards: "--discard-untrimmed --rename='{id}_{r1.adapter_name} {comment}'"
    log:
        os.path.join(config['dir']['log'], "{sample}_cutadapt_5end.log")
    threads:
        config['threads']
    wrapper:
        "0.79.0/bio/cutadapt/pe"


rule cutadapt_5end_rename_R2:
    input:
        os.path.join(config['dir']['raw'], '{sample}_1.fq.gz'),
        os.path.join(config['dir']['raw'], '{sample}_2.fq.gz')
    output:
        fastq1 = os.path.join(config['dir']['results']['fq'], '{sample}_5end_1.fq.gz'),
        fastq2 = os.path.join(config['dir']['results']['fq'], '{sample}_5end_2.fq.gz'),
        qc = os.path.join(config['dir']['qc'], "{sample}_cutadapt_5end.txt")
    params:
        adapters = "-g file:{R1} -G file:{R2}".format(R1=config['adapters']['R1']['5end'],
                                                      R2=config['adapters']['R2']['5end']),
        extra = lambda wildcards: "--discard-untrimmed --rename='{id}_{r2.adapter_name} {comment}'"
    log:
        os.path.join(config['dir']['log'], "{sample}_cutadapt_5end.log")
    threads:
        config['threads']
    wrapper:
        "0.79.0/bio/cutadapt/pe"


rule cutadapt_3end_se_R1:
    input:
        os.path.join(config['dir']['results']['fq'], '{sample}_5end_1.fq.gz')
    output:
        fastq=os.path.join(config['dir']['results']['fq'], '{sample}_3end_se_1.fq.gz'),
        qc=os.path.join(config['dir']['qc'], "{sample}_cutadapt_3end_se_R1.txt")
    params:
        adapters="-a file:{R1}".format(R1=config['adapters']['R1']['3end_se']),
        extra="",
    log:
        os.path.join(config['dir']['log'], "{sample}_cutadapt_3end.log")
    threads: config['threads']
    wrapper:
        "0.79.0/bio/cutadapt/se"


rule cutadapt_3end:
    input:
        os.path.join(config['dir']['results']['fq'], '{sample}_5end_1.fq.gz'),
        os.path.join(config['dir']['results']['fq'], '{sample}_5end_2.fq.gz')
    output:
        fastq1 = os.path.join(config['dir']['results']['fq'], '{sample}_3end_1.fq.gz'),
        fastq2 = os.path.join(config['dir']['results']['fq'], '{sample}_3end_2.fq.gz'),
        qc = os.path.join(config['dir']['qc'], "{sample}_cutadapt_3end.txt")
    params:
        adapters = "-a file:{R1} -A file:{R2}".format(R1=config['adapters']['R1']['3end'],
                                                      R2=config['adapters']['R2']['3end']),
        extra = "-m 30",
    log:
        os.path.join(config['dir']['log'], "{sample}_cutadapt_3end.log")
    threads: config['threads']
    wrapper:
        "0.79.0/bio/cutadapt/pe"


rule cutadapt_separate_by_R2:
    input:
        os.path.join(config['dir']['results']['fq'], '{sample}_3end_1.fq.gz'),
        os.path.join(config['dir']['results']['fq'], '{sample}_3end_2.fq.gz'),
    output:
        fastq1 = os.path.join(config['dir']['results']['fq'], '{sample}_spikeins_1.fq.gz'),
        fastq2 = os.path.join(config['dir']['results']['fq'], '{sample}_spikeins_2.fq.gz'),
        fastq3 = os.path.join(config['dir']['results']['fq'], '{sample}_genome_1.fq.gz'), # fake output
        fastq4 = os.path.join(config['dir']['results']['fq'], '{sample}_genome_2.fq.gz'), # fake output
        qc = os.path.join(config['dir']['qc'], "{sample}_cutadapt_separate.txt")
    params:
        adapters = "-G file:{sep}".format(sep=config['adapter_spikeins']),
        extra = lambda wildcards: "-e 0 -o 8 --action=none --untrimmed-o {0[0]} --untrimmed-p {0[1]}".format(
                        expand(os.path.join(config['dir']['results']['fq'],
                                            "{sample}_genome_{num}.fq.gz"),
                        sample=wildcards.sample, num=[1, 2])
        )
    log:
        os.path.join(config['dir']['log'], "{sample}_cutadapt_separate.log")
    threads:
        config['threads']
    wrapper:
        "0.79.0/bio/cutadapt/pe"


rule cutadapt_separate_by_R1:
    input:
        os.path.join(config['dir']['results']['fq'], '{sample}_3end_1.fq.gz'),
        os.path.join(config['dir']['results']['fq'], '{sample}_3end_2.fq.gz'),
    output:
        fastq1 = os.path.join(config['dir']['results']['fq'], '{sample}_spikeins_1.fq.gz'),
        fastq2 = os.path.join(config['dir']['results']['fq'], '{sample}_spikeins_2.fq.gz'),
        fastq3 = os.path.join(config['dir']['results']['fq'], '{sample}_genome_1.fq.gz'), # fake output
        fastq4 = os.path.join(config['dir']['results']['fq'], '{sample}_genome_2.fq.gz'), # fake output
        qc = os.path.join(config['dir']['qc'], "{sample}_cutadapt_separate.txt")
    params:
        adapters = "-g file:{sep}".format(sep=config['adapter_spikeins']),
        extra = lambda wildcards: "-e 0 -o 8 --action=none --untrimmed-o {0[0]} --untrimmed-p {0[1]}".format(
                        expand(os.path.join(config['dir']['results']['fq'],
                                            "{sample}_genome_{num}.fq.gz"),
                        sample=wildcards.sample, num=[1, 2])
        )
    log:
        os.path.join(config['dir']['log'], "{sample}_cutadapt_separate.log")
    threads:
        config['threads']
    wrapper:
        "0.79.0/bio/cutadapt/pe"


rule bowtie2:
    input:
        sample = define_bowtie2_input(config['adapter_spikeins'])
    output:
        os.path.join(config['dir']['results']['bam'], "{sample}_{tag}_aligned.bam")
    log:
        os.path.join(config['dir']['qc'], "{sample}_{tag}_aligned.txt")
    params:
        index = lambda wildcards: config['refs'][wildcards.tag]['index']['bowtie2'],
        extra = "-X 2000 --no-discordant --no-mixed --no-unal"
    threads: config['threads']  # Use at least two threads
    wrapper:
        "0.79.0/bio/bowtie2/align"


rule filter_by_mapQ:
    input:
        os.path.join(config['dir']['results']['bam'], "{sample}_genome_aligned.bam")
    output:
        os.path.join(config['dir']['results']['bam'], "{sample}_genome_filtered.bam")
    params:
        extra = "-Shb -f 0x2 -q {mapQ}".format(mapQ=config['params']['mapQ'])
    log:
        os.path.join(config['dir']['log'], "{sample}_genome_samtools_view.log")
    threads: config['threads']
    wrapper:
        "0.79.0/bio/samtools/view"


rule filter_by_mapQ_and_chromosome:
    input:
        os.path.join(config['dir']['results']['bam'], "{sample}_genome_aligned.bam")
    output:
        os.path.join(config['dir']['results']['bam'], "{sample}_genome_filtered.bam")
    params:
        extra="-Shb -f 0x2 -q {mapQ} -L {chroms}".format(mapQ=config['params']['mapQ'], chroms=config['refs']['genome']['chromosome'])
    log:
        os.path.join(config['dir']['log'], "{sample}_genome_samtools_view.log")
    threads: config['threads']
    wrapper:
        "0.79.0/bio/samtools/view"


rule samtools_nsort:
    input:
        lambda wildcards: os.path.join(config['dir']['results']['bam'], "{a}_{b}.bam").format(
            a=wildcards.sample,
            b=f"{wildcards.tag}_aligned" if wildcards.tag == 'spikeins' else f"{wildcards.tag}_filtered")
    output:
        os.path.join(config['dir']['results']['bam'], "{sample}_{tag}_nsorted.bam")
    params:
        extra = "-m 4G -n",
        tmp_dir = "/tmp/"
    threads:
        config['threads']     # This value - 1 will be sent to -@.
    log:
        os.path.join(config['dir']['log'], "{sample}_{tag}_samtools_nsort.log")
    wrapper:
        "0.79.0/bio/samtools/sort"


rule bedtools_bamtobed:
    input:
        os.path.join(config['dir']['results']['bam'], "{sample}_{tag}_nsorted.bam")
    output:
        os.path.join(config['dir']['results']['bedpe'], "{sample}_{tag}.bedpe")
    params:
        "-bedpe -mate1"
    log:
        os.path.join(config['dir']['log'], "{sample}_{tag}_bedtools_bamtobed.log")
    conda:
        "env/bedtools.yaml"
    shell:
        "bedtools bamtobed {params} -i {input} | bedtools sort -i stdin > {output} 2> {log}"


rule count_fragments:
    input:
        bedpe = os.path.join(config['dir']['results']['bedpe'], "{sample}_{tag}.bedpe")
    output:
        counts = os.path.join(config['dir']['results']['bed'], '{sample}_{tag}.bed'),
    log:
        os.path.join(config['dir']['log'], "{sample}_{tag}_count_fragments.log")
    params:
        sample = lambda wildcards: wildcards.sample,
        reads_renamed = config['adapters']['reads_renamed'],
        strand_as_R1 = config['adapters']['strand_as_R1'],
    conda:
        "env/R.yaml"
    threads: 8
    script:
        "scripts/count_fragments.R"


rule multiqc:
    input:
        #define_final_output(config, subset, config['dir']['qc'], '_{num}_fastqc.zip'),
        #define_final_output(config, subset, config['dir']['qc'], '_cutadapt_5end.txt'),
        #define_final_output(config, subset, config['dir']['qc'], '_cutadapt_3end.txt'),
        #define_final_output(config, subset, config['dir']['qc'], '_cutadapt_separate.txt'),
        #define_final_output(config, subset, config['dir']['qc'], '_{tag}_aligned.txt'),
    output:
        #os.path.join(config['dir']['qc'], "multiqc.html")
    log:
        #os.path.join(config['dir']['log'], "multiqc.log")
    wrapper:
        "0.79.0/bio/multiqc"
