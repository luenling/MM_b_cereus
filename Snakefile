import glob as glob
import re as re
# specify name of external config-file
configfile: "MM_b_cereus/extraconfs.yaml"

# import custom global parameter from config-file
DATASETS_RNA = config["datasets_rna"]
DATASETS_DNA = config["datasets_dna"]
DATASETS = DATASETS_RNA + DATASETS_DNA
CONDITIONS = config["conditions"]
THREADS = config["threads"]
#dgeaCmpsString = config["dgea-cmps"]
#dgeaCmpsList = dgeaCmpsString.split("::")
# readtype can be single or paired-end
readType = config["readType"]

# list of rules which are not deployed to slurm
localrules: all, fastqc, quant, star, salmon_index, coco, clean


# final target rule to produce all sub targets
rule all:
    input:
        fastqc = "logs/fastqc.done",
        #quant  = "logs/quant.done",
        #star   = "logs/star.done",
        #coco   = "logs/coco.done"
#        enrich = "logs/enrichmentBrowser.done",
#        master = "logs/master.done"


# subtarget: bam to fastq conversion
rule fastqc:
    input:
        fastqtrim = expand("qc/fastqc_trimmed/{dataset}_fastqc.html", dataset=DATASETS),
        fastqraw = expand("qc/fastqc_raw/{dataset}_fastqc.html", dataset=DATASETS),
    output: "logs/fastqc.done"
    shell:
        '''
        touch {output}
        '''

# subtarget: read per transcript quantification
rule quant:
    input:
        bam = expand("quants/{dataset}/alignments.bam", dataset=DATASETS)
    output: "logs/quant.done"
    shell:
        '''
        touch {output}
        '''

rule star:
    input:
        bam = expand("STAR/{dataset}/{dataset}_Aligned.sortedByCoord.out.bam", dataset=DATASETS)
    output: "logs/star.done"
    shell:
        '''
        touch {output}
        '''
rule picardqc:
    input:
        picardqc = expand("qc/picard/{dataset}.RNASeqMetrics.txt", dataset=DATASETS)
    output: "logs/picardqc.done"
    shell:
        '''
        touch {output}
        '''

###
######
###

def get_fastq_fn(wildcards):
    fastqfn = glob.glob(f"raw/*{wildcards.sample}*.fastq.gz")
    return(fastqfn)

# step 1: check read quality of raw reads
rule fastqc_check:
    input:
        read1 = get_bam_fn
    output:
        fastqc = "qc/fastqc_raw/{sample}_fastqc.html"
    threads: THREADS
    params:
        dir = "qc/fastqc_raw",
        sample = "{sample}"
    # run:
    #     print(f"{params.sample}")
    #     print(glob.glob(f"raw/*{params.sample}*.bam"))
    #     print(f"{input.reads}")
    shell:
        '''
        fastqc -Xmx4g --outdir {params.dir} --thread {threads} {input.reads}  > /dev/null 2> /dev/null
        '''


# step 2: clip/trimm reads with bbduk
rule bbduk:
    input:
        reads1 = "raw/{sample}_L001_R1_001.fastq.gz"
        reads2 = "raw/{sample}_L001_R2_001.fastq.gz"
    output:
        fastq1 = "fastq/{sample}_1.fastq.gz"
        fastq2 = "fastq/{sample}_2.fastq.gz"
    log:
        err = "logs/bbduk/{sample}.log"
    threads: THREADS
    params:
        bbduk = config['bbduk'],
        bbduk_parms = config['bbduk-params'],
        bbduk_adapt = config['bbduk-adapt']
    shell:
        "{params.bbduk} -Xmx6g in1={input.reads1} in2={input.reads2} "
        "out1={output.fastq1} out2={output.fastq2} " 
        "ref={params.bbduk_adapt} {params.bbduk_parms} > /dev/null 2> {log.err}"

# step 3: check read quality of trimmed reads
rule fastqc_recheck:
    input:
        fastq = "fastq/{sample}.fastq.gz"
    output:
        fastqc = "qc/fastqc_trimmed/{sample}_fastqc.html"
    threads: THREADS
    shell:
        '''
        fastqc  -Xmx6g --outdir qc/fastqc_trimmed/ --thread {threads} {input.fastq}  > /dev/null 2> /dev/null
        '''


# step 6: hard link prepared salmon index (v1.1.0)
rule salmon_index:
    input:
        idx = config["salmon-idx"],
        t2g = config["salmon-t2g"]
    output:
        idx = directory("auxData/salmon.idx"),
        t2g = "auxData/t2g.csv"
    threads: THREADS
    shell:
        '''
        mkdir -p auxData/
        cp -r -l {input.idx} {output.idx}
        cp -r -l {input.t2g} {output.t2g}
        '''

# step 7: quantify tanscript abundance with salmon (v1.1.0)
rule salmon_quant:
    input:
        mate = "fastq/{sample}.fastq.gz",
        idx   = "auxData/salmon.idx",
        t2g   = "auxData/t2g.csv"
    output:
        sf  = "quants/{sample}/quant.sf",
        gsf = "quants/{sample}/quant.genes.sf",
        sam = temp("quants/{sample}/alignments.sam"),
    log:
        err = "logs/salmon_quant/{sample}.err",
        out = "logs/salmon_quant/{sample}.out"
    threads: 6
    resources:
        time_long = "04:00:00",
        mem_mb = 64000
    params:
        dir = "quants/{sample}",
        params = config["salmon-params"],
        libtype = config["salmon-libtype"]
    shell:
        "module load Salmon/1.4.0-gompi-2020b\n"
        "salmon quant --quiet --threads {threads} {params.params} "
        "-i {input.idx} --libType {params.libtype} -r {input.mate} "
        "--geneMap {input.t2g} --writeMappings={output.sam} "
        "--output {params.dir}  > {log.out} 2> {log.err}"

# step 8: generate sorted bam file from salmon generated sam file
rule sam2sortedbam_salmon:
    input: "quants/{sample}/alignments.sam"
    output:
        bam = "quants/{sample}/alignments.bam",
        bai = "quants/{sample}/alignments.bam.bai"
    log: "logs/sam2sortedbam/{sample}.log"
    threads: 4
    shell:
        '''
        samtools sort -m 2G --threads {threads} -o {output.bam} {input} 2> {log}
        samtools index {output.bam} 2>> {log}
        '''
# step 9: map with star

rule star_map:
    input:
        fastq = "fastq/{sample}.fastq.gz"
    output:
        bam  = "STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
    log: "logs/STAR/{sample}.log"
    threads: 15
    resources:
        time_long = "04:00:00",
        mem_mb = 64000
    params:
        outprefix = "STAR/{sample}/{sample}_",
        idx  = directory(config["mm10_dir"] + "/" + config["mm10_star_idx"])
    shell:
        "module load STAR/2.7.9a-GCC-10.3.0\n"
        "STAR --runThreadN {threads} "
        " --outFileNamePrefix {params.outprefix} "
        "--genomeDir {params.idx} "
        "--readFilesIn {input.fastq} "
        "--readFilesCommand  gunzip -c "
        "--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts > {log} 2>&1\n"
        "samtools index {output} 2>> {log}"

## not used, jut for recreating the index
rule star_index:
    input:
        fasta = config["mm10_fasta"],
        gtf   = config["mm10_gtf"]
    output:
        index_dir = directory(config["mm10_dir"] + "/" + config["mm10_star_idx"])
    log: "logs/STAR/index.log"
    threads: 25
    resources:
        time_long = "04:00:00",
        mem_mb = 64000
    params:
        genomeDir = config["mm10_dir"],
        idxDir    = config["mm10_star_idx"]
    shell:
        "module load STAR/2.7.9a-GCC-10.3.0\n"
        "STAR --runThreadN {threads} --runMode genomeGenerate "
        "--genomeDir {params.genomeDir} "
        "--genomeFastaFiles {input.fasta} "
        "--sjdbGTFfile {input.gtf} "
        "--sjdbOverhang 99 2>&1 > {log}"


rule picard_RNAqc:
    input: "STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
    output: "qc/picard/{sample}.RNASeqMetrics.txt"
    params:
        picard = config["picard"],
        refflat =  config["refflat"],
        ribos = config["ribos"]
    log: "logs/picard/{sample}.log"
    threads: THREADS
    shell:
        "java -Xmx6g -jar {params.picard} CollectRnaSeqMetrics "
        "I={input} O={output} "
        "REF_FLAT={params.refflat} "
        "STRAND=NONE "
        "RIBOSOMAL_INTERVALS={params.ribos} "
        "VALIDATION_STRINGENCY=LENIENT VERBOSITY=WARNING > {log} 2>&1 "

## get samtools stats for qc
rule samtool_stats:
    input : "STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
    output : "qc/samtools/{sample}_sam_stats.txt"
    params :
        fasta = config["mm10_fasta"]
    threads : 3
    log : "logs/samtools/{sample}.log"
    shell:
        "samtools stats -r {params.fasta} --threads {threads} "
        "{input} > {output} 2> {log}"

rule multiqc:
    input : expand("qc/picard/{dataset}.RNASeqMetrics.txt", dataset = DATASETS),
            expand("qc/samtools/{dataset}_sam_stats.txt",  dataset = DATASETS)
    output : "multiqc/multiqc_report.html"
    threads : THREADS
    log : "logs/multiqc.log"
    shell:
        "module load MultiQC\n"
        "multiqc --no-ansi -f -o multiqc qc logs STAR quants FeatureCount > {log} 2>&1"

## feature count just included for reverse stranded single end
rule feature_count:
    input : expand("STAR/{dataset}/{dataset}_Aligned.sortedByCoord.out.bam", dataset=DATASETS)
    output : "FeatureCount/count_table.tsv"
    params :
        gtf = config["mm10_gtf"]
    threads : 10
    log : "logs/featurecounts.log"
    shell:
        "module load Subread\n"
        "featureCounts -T {threads} -d 30 -s 2 -Q 1 "
        "-g gene_id -F GTF -a {params.gtf} "
        "--minOverlap 10 --largestOverlap "
        "-o {output} {input} 2> {log}"

rule coco:
    input: expand("Coco/{dataset}_counts.tsv", dataset=DATASETS)
    output: "logs/coco.done"
    shell:
        '''
        touch {output}
        '''

## feature count just included for reverse stranded single end
rule coco_count:
    input : "STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
    output : "Coco/{sample}_counts.tsv"
    params :
        coco = config["coco"],
        gtf  =  config["coco_gtf"]
    threads : 3
    log : "logs/coco/coco_{sample}.log"
    shell:
        "module load Subread\n"
        "{params.coco} correct_count -s 2 -t {threads} "
        " {params.gtf} {input} {output} 2> {log}"


###
######
###

# rule to remove all workflow intermediate and results files, but the logs
rule clean:
    threads: THREADS
    shell:
        '''
        rm -rf STAR/
        rm -rf fastq/
        rm -rf quants/
        rm -rf qc/
        rm -rf dgea/
        rm -rf auxData/
        rm -rf enrichmentBrowser/
        rm -rf masterTable/
        '''


## mmquant does some weird clustering of ambiguously mapped reads  not that useful after all, coco seems to be better in handling multimapped reads
rule mmquant:
    input : expand("STAR/{dataset}/{dataset}_Aligned.sortedByCoord.out.bam", dataset=DATASETS)
    output :
        countfile = "mmquant/mmquant_count_table.tsv",
        statfile =  "mmquant/mmquant_stats.txt"
    params :
        gtf = config["mm10_gtf"],
        mmquant = config["mmquant"]
    threads : 15
    log : "logs/mmquant.log"
    shell:
        "{params.mmquant} -t {threads} -f bam -s R -d 15 -l 10 "
        "-O {output.statfile} -a {params.gtf} "
        "-o {output.countfile} -r {input} 2> {log}"




## not yet working - might not be worth it, all bed files need to be unzipped
rule rseqc:
    input: "STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        gbcov =  "qc/rseqc/{sample}/{sample}.geneBodyCoverage.txt",
        juncanno = "qc/rseqc/{sample}/{sample}.junction_annotation.txt",
        readdist = "qc/rseqc/{sample}/{sample}.read_distribution.txt",
        infexp = "qc/rseqc/{sample}/{sample}.infer_experiment.txt",
        tin = "qc/rseqc/{sample}.tin.txt"
    params:
        hk = config["rseqc-hkgenes"],
        ribo =  config["rseqc-ribogenes"],
        genes =  config["rseqc-genes"],
    threads: 2
    log : "logs/rseqc/{sample}.log"
    shell:
        '''
        module load RSeQC
        tin.py -i {input} -r {params.genes} > {output.infexp} 2>> {log}
        junction_annotation.py -i {input} -o qc/rseqc/{wildcards.sample}/{wildcards.sample} -r {params.genes} 2> {output.juncanno}
        geneBody_coverage.py -i {input} -o qc/rseqc/{wildcards.sample}/{wildcards.sample} -r {params.genes} 2> qc/rseqc/{wildcards.sample}/{wildcards.sample}.gbcoverage_error.txt
        '''


# step 9: make differential gene experession analysis
# run external R-script (scripts/DESeq2.R)
# see script for requirments concerning file name pattern
rule deseq2_dgea:
    input:
        quants = expand("quants/{dataset}/quant.sf", dataset=DATASETS),
        g2t = "auxData/t2g.csv",
        bin = "scripts/DESeq2_multiComp.R"
    output:
        pca  = "dgea/DESeq2_pca.pdf",
        heat = "dgea/DESeq2_heatmap.pdf",
        csv  = expand("dgea/DESeq2_{cmp}.csv", cmp = dgeaCmpsList),
        cnts = "dgea/DESeq2_counts.csv",
        tpm  = expand("dgea/DESeq2_TPM_{cond}.csv", cond = CONDITIONS),
        cnt  = expand("dgea/DESeq2_counts_{cmp}.csv", cmp = dgeaCmpsList),
        ma   = expand("dgea/DESeq2_MAplot_{cmp}.pdf", cmp = dgeaCmpsList),
        si   = "dgea/RsessionInfo.txt"
    threads: THREADS
    params:
        cmps = dgeaCmpsString
    log:
        err = "logs/deseq2_dgea/DESeq2.err",
        out = "logs/deseq2_dgea/DESeq2.out"
    shell:
        '''
        module purge
        module load R/4.2.1-foss-2022a
        Rscript --no-save --no-restore {input.bin} {input.g2t} {params.cmps} {input.quants} > {log.out} 2> {log.err}
        '''

# step 12: make functional enrichment analysis with enrichmentBrowser
## extract read counts from salmon
rule enrichmentPrep:
    input:
        quants = expand("dgea/DESeq2_counts_{cmp}.csv", cmp = dgeaCmpsList),
        bin = "scripts/enrichmentBrowser_multiComp.pl"
    output:
        cnt = expand("enrichmentBrowser/data/cnts_{cmp}.tab", cmp = dgeaCmpsList),
        row = expand("enrichmentBrowser/data/rows_{cmp}.tab", cmp = dgeaCmpsList),
        col = expand("enrichmentBrowser/data/cols_{cmp}.tab", cmp = dgeaCmpsList)
    threads: THREADS
    params:
        cmps = dgeaCmpsString
    log:
        err = "logs/enrichmentBrowser/enrichmentPrep.err"
    shell:
        '''
        perl {input.bin} {input.quants} 2> {log.err}
        '''

# step 14: run enrichmentBrowser script
rule enrichmentRun:
    input:
        dat  = "enrichmentBrowser/data/cnts_{cmp}.tab",
        row  = "enrichmentBrowser/data/rows_{cmp}.tab",
        col  = "enrichmentBrowser/data/cols_{cmp}.tab",
        stat = "dgea/DESeq2_{cmp}.csv",
        bin  = "scripts/enrichmentBrowser_multiComp.R",
    output:
        res = "enrichmentBrowser/{cmp}/go/index.html"
    params:
        pval = config["enrich-pval"],
        lfc  = config["enrich-lfc"],
        perm = config["enrich-perm"],
        org  = config["enrich-org"],
        ids  = config["enrich-ids"],
        anno = config["enrich-anno"],
        show = config["enrich-show"],
        gsm  = config["enrich-gsMethod"],
        grnm = config["enrich-grnMethod"]
    log:
        out = "logs/enrichmentBrowser/enrichmentRun_{cmp}.log",
        err = "logs/enrichmentBrowser/enrichmentRun_{cmp}.err"
    threads: THREADS
    shell:
        '''
        module load conda
        # conda create -n enrichmentbrowser -c bioconda bioconductor-enrichmentbrowser
        conda activate enrichmentbrowser
        Rscript --no-save --no-restore {input.bin} {params.pval} {params.lfc} {params.perm} {params.org} {params.ids} {params.anno} {params.show} {params.gsm} {params.grnm} {input.stat}  > {log.out} 2> {log.err}
        '''

# step 15: copy auxData needed to comile master table
rule master_aux:
    input:
        biomart = config["masterTab-idMap"],
    output:
        biomart = "auxData/biomart.tsv.gz",
    threads: THREADS
    shell:
        '''
        cp {input.biomart} {output.biomart}
        '''
