configfile: "config.yaml"

SAMPLES, = glob_wildcards(config["DIRECTORY_FQ"] + "{sample}.fastq")

rule all:
    input:
        "mappen.aangemaakt",
        "SNP.treefile",
        "flye_tree.newick"


# Benodigde mappen aanmaken
rule create_dirs:
    output:
        "mappen.aangemaakt"
    shell:
        """
        mkdir -p qc_results qc_results/old qc_results/edited edited_fastq BAM pileup bcf consensus flye hashing tmp benchmarks
        touch mappen.aangemaakt
        """


# Deelopdracht 1
rule quality_control_old:
    input:
        fastq=config["DIRECTORY_FQ"] + "{sample}.fastq"
    output:
        json="qc_results/old/{sample}_qc.json"
    benchmark:
        "benchmarks/qc_{sample}.txt"
    script:
        "quality_control.py"

# Deelopdracht 2
rule filter_trimming:
    input:
        fastq=config["DIRECTORY_FQ"] + "{sample}.fastq",
        qc="qc_results/old/{sample}_qc.json"
    output:
        fastq_edited="edited_fastq/{sample}.fastq"
    benchmark:
        "benchmarks/trimming_{sample}.txt"
    script:
        "filter_trimming.py"

# Deelopdracht 3
rule quality_control_edited:
    input:
        fastq="edited_fastq/{sample}.fastq"
    output:
        json="qc_results/edited/{sample}_qc.json"
    script:
        "quality_control.py"

# Deelopdracht 4
rule align_minimap2:
    input:
        fastq="edited_fastq/{sample}.fastq",
        json="qc_results/edited/{sample}_qc.json",
        refgenome=config["REFGENOME"]
    output:
        sam=temp("BAM/{sample}.sam")
    threads: 4
    shell:
        """
        minimap2 -t {threads} -x map-ont -a {input.refgenome} {input.fastq} > {output.sam}
        """

rule sam_to_bam:
    input:
        sam="BAM/{sample}.sam"
    output:
        bam=temp("BAM/{sample}.bam"),
        sorted_bam="BAM/{sample}.sorted.bam"
    threads: 4
    shell:
        """
        samtools view -@ {threads} -bS {input.sam} > {output.bam}
        samtools sort -@ {threads} -o {output.sorted_bam} {output.bam}
        samtools index {output.sorted_bam}
        """

# Deelopdracht 5
rule faidx:
    input:
        refgenome=config["REFGENOME"]
    output:
        refgenome_fai=config["REFGENOME_BASE_NAME"] + ".fai",
        refgenome_base=config["REFGENOME_BASE_NAME"]
    shell:
        """
        cp {input.refgenome} .
        samtools faidx {output.refgenome_base}
        """

rule mpileup:
    input:
        bam="BAM/{sample}.sorted.bam",
        refgenome=config["REFGENOME_BASE_NAME"],
        refgenome_fai=config["REFGENOME_BASE_NAME"] + ".fai"
    output:
        pileup="pileup/{sample}.pileup"
    shell:
        "bcftools mpileup -f {input.refgenome} -o {output.pileup} {input.bam}"

rule call_variants:
    input:
        pileup="pileup/{sample}.pileup"
    output:
        vcf="bcf/{sample}.vcf"
    shell:
        "bcftools call -mv -Ov -P 0.9 -o {output.vcf} {input.pileup} "

rule filter_snps:
    input:
        vcf="bcf/{sample}.vcf"
    output:
        snps_vcf="bcf/{sample}_snps.vcf"
    shell:
        "bcftools filter -i 'TYPE=\"snp\"' {input.vcf} -o {output.snps_vcf}"

rule compress_to_bcf:
    input:
        snps_vcf="bcf/{sample}_snps.vcf"
    output:
        bcf="bcf/{sample}.bcf"
    shell:
        "bcftools view -Ob -o {output.bcf} {input.snps_vcf}"

rule index_bcf:
    input:
        bcf="bcf/{sample}.bcf"
    output:
        bcf_index="bcf/{sample}.bcf.csi"
    shell:
        """
        bcftools index {input.bcf}
        """

rule consensus:
    input:
        bcf="bcf/{sample}.bcf",
        bcf_index="bcf/{sample}.bcf.csi",
        refgenome=config["REFGENOME"],
        refgenome_fai=config["REFGENOME_BASE_NAME"] + ".fai"
    output:
        fasta="consensus/{sample}_consensus.fasta"
    shell:
        """
        bcftools consensus -f {input.refgenome} {input.bcf} -o {output.fasta}
        """

# Deelopdracht 6
rule modify_headers:
    input:
        fasta="consensus/{sample}_consensus.fasta"
    output:
        modified_fasta="consensus/{sample}_modified.fasta"
    run:
        with open(input.fasta) as infile, open(output.modified_fasta, 'w') as outfile:
            for i, line in enumerate(infile):
                if line.startswith('>'):
                    line = f">{wildcards.sample}_{i // 2}\n"
                outfile.write(line)

rule combine_fastas:
    input:
        expand("consensus/{sample}_modified.fasta", sample=SAMPLES)
    output:
        combined_fasta="consensus/combined.fasta"
    run:
        with open(output.combined_fasta, 'w') as outfile:
            for fasta in input:
                with open(fasta) as infile:
                    outfile.write(infile.read())

rule run_iqtree:
    input:
        aligned_fasta="consensus/combined.fasta"
    output:
        tree_file="SNP.treefile"
    shell:
        """
        iqtree2 -s {input} -nt AUTO -bb 1000 -alrt 1000 -pre SNP
        """


# Deelopdracht 7
rule run_flye:
    input:
        fastq="edited_fastq/{sample}.fastq"
    output:
        directory("flye/{sample}/")
    resources:
        lock=1
    shell:
        "flye --nano-raw {input.fastq} --out-dir {output} --threads 4"


# Deelopdracht 8
rule hashing:
    input:
        assembly="flye/{sample}"
    output:
        hash_sketch="hashing/{sample}.sig"
    shell:
        "sourmash sketch dna -p k=21,scaled=100 --output {output.hash_sketch} {input.assembly}/assembly.fasta"

rule distance_matrix:
    input:
        expand("hashing/{sample}.sig", sample=SAMPLES)
    output:
        comparison="comparison.csv"
    shell:
        "sourmash compare hashing/*.sig --csv {output.comparison}"

rule make_flye_tree:
    input:
        comparison="comparison.csv"
    output:
        tree="flye_tree.newick"
    shell:
        "python " + config["DIST_TO_NEWICK"] + " {input.comparison} {output.tree}"
